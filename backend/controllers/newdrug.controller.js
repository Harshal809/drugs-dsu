// backend/controllers/finalnewdrug.controller.js
import { spawn } from 'child_process';
import path from 'path';
import { fileURLToPath } from 'url';
import os from 'os';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// Resolve repo root from backend/controllers -> repo root
const repoRoot = path.resolve(__dirname, '..', '..');

// Path to your Python MCP server entry
const pyEntry = path.join(repoRoot, 'mcpserver', 'main.py');

// Choose python binary (allow override through env)
const PYTHON_BIN = process.env.PYTHON_BIN || (os.platform() === 'win32' ? 'python' : 'python3');

/***********************
 * Minimal MCP Client
 ***********************/
class MCPClient {
  constructor() {
    this.proc = null;
    this.nextId = 1;
    this.pending = new Map(); // id -> {resolve, reject, timeout}
    this.buffer = Buffer.alloc(0);
    this.ready = null;
    this._starting = false;
  }

  async start() {
    if (this.ready) return this.ready;
    if (this._starting) return this._starting;

    this._starting = new Promise((resolve, reject) => {
      try {
        this.proc = spawn(PYTHON_BIN, [pyEntry], {
          cwd: repoRoot,
          stdio: ['pipe', 'pipe', 'pipe'],
          env: { ...process.env },
        });

        this.proc.on('error', (err) => {
          reject(new Error(`Failed to start MCP server: ${err.message}`));
        });

        this.proc.stderr.on('data', (data) => {
          // Optional: log for debugging
          // console.error('[MCP STDERR]', data.toString());
        });

        this.proc.stdout.on('data', (chunk) => {
          this._onData(chunk);
        });

        this.proc.on('exit', (code, signal) => {
          const err = new Error(`MCP server exited with code=${code} signal=${signal || ''}`);
          for (const [, p] of this.pending.entries()) {
            clearTimeout(p.timeout);
            p.reject(err);
          }
          this.pending.clear();
          this.proc = null;
          this.buffer = Buffer.alloc(0);
          this.ready = null;
          this._starting = false;
        });

        // Once process is up, send initialize
        this.initialize()
          .then(() => {
            this.ready = Promise.resolve(true);
            resolve(true);
          })
          .catch(reject);
      } catch (e) {
        reject(e);
      }
    });

    return this._starting;
  }

  stop() {
    if (this.proc) {
      this.proc.kill();
      this.proc = null;
    }
  }

  _onData(chunk) {
    this.buffer = Buffer.concat([this.buffer, chunk]);

    // Parse LSP-style headers
    while (true) {
      const headerEnd = this.buffer.indexOf('\r\n\r\n');
      if (headerEnd === -1) break;

      const headerStr = this.buffer.slice(0, headerEnd).toString('utf8');
      const headers = {};
      for (const line of headerStr.split('\r\n')) {
        const idx = line.indexOf(':');
        if (idx > -1) {
          const k = line.slice(0, idx).trim().toLowerCase();
          const v = line.slice(idx + 1).trim();
          headers[k] = v;
        }
      }

      const contentLength = parseInt(headers['content-length'] || '0', 10);
      const totalLen = headerEnd + 4 + contentLength;
      if (this.buffer.length < totalLen) break;

      const bodyBuf = this.buffer.slice(headerEnd + 4, totalLen);
      this.buffer = this.buffer.slice(totalLen);

      let msg = null;
      try {
        msg = JSON.parse(bodyBuf.toString('utf8'));
      } catch {
        continue;
      }

      this._dispatch(msg);
    }
  }

  _dispatch(msg) {
    const id = msg.id;
    if (id == null) {
      // unsolicited or error without id
      return;
    }
    const pending = this.pending.get(id);
    if (!pending) return;

    clearTimeout(pending.timeout);
    this.pending.delete(id);

    if (msg.error) {
      pending.reject(new Error(`${msg.error.code || ''} ${msg.error.message || 'Unknown MCP error'}`));
    } else {
      pending.resolve(msg.result);
    }
  }

  _send(payload, timeoutMs = 60000) {
    return new Promise((resolve, reject) => {
      if (!this.proc || !this.proc.stdin.writable) {
        return reject(new Error('MCP server is not running'));
      }

      const id = payload.id ?? this.nextId++;
      payload.id = id;

      const data = Buffer.from(JSON.stringify(payload), 'utf8');
      const header = Buffer.from(`Content-Length: ${data.length}\r\n\r\n`, 'utf8');
      const toWrite = Buffer.concat([header, data]);

      const to = setTimeout(() => {
        this.pending.delete(id);
        reject(new Error(`MCP request timeout (id=${id}, method=${payload.method})`));
      }, timeoutMs);

      this.pending.set(id, { resolve, reject, timeout: to });

      try {
        this.proc.stdin.write(toWrite);
      } catch (e) {
        clearTimeout(to);
        this.pending.delete(id);
        reject(e);
      }
    });
  }

  async initialize() {
    await this.start();
    return this._send(
      {
        method: 'initialize',
        params: {},
      },
      15000
    );
  }

  async toolsCall(name, args = {}, timeoutMs = 90000) {
    await this.start();
    const res = await this._send(
      {
        method: 'tools/call',
        params: { name, arguments: args },
      },
      timeoutMs
    );
    return res; // Typically { content: [...] }
  }


  async analyzeTargets(symptoms, max_targets = 3) {
    return this.toolsCall('analyze_targets', {
      symptoms: String(symptoms || ''),
      max_targets: Number(max_targets || 3),
    });
  }

  async ligandsForProtein({ protein_id = '', protein_name = '', max_ligands = 5, include_admet = true }) {
    return this.toolsCall('get_ligands_for_protein', {
      protein_id: String(protein_id || ''),
      protein_name: String(protein_name || ''),
      max_ligands: Number(max_ligands || 5),
      include_admet: Boolean(include_admet),
    });
  }

  async health() {
    return this.toolsCall('health', {});
  }
}

export const mcp = new MCPClient();

// Graceful stop on process exit
process.on('exit', () => mcp.stop());
process.on('SIGINT', () => { mcp.stop(); process.exit(0); });
process.on('SIGTERM', () => { mcp.stop(); process.exit(0); });

// Optional: pre-warm so first request is fast (non-blocking)
mcp.start().catch((e) => console.error('MCP start failed:', e));

/***********************
 * Controllers
 ***********************/

// Helper to unwrap MCP tool JSON content
function extractToolJson(result) {
  try {
    const item = (result?.content || []).find((c) => c?.type === 'json');
    return item?.json ?? null;
  } catch {
    return null;
  }
}

// POST /predictDisease/:id
// Body: { symptoms: "fever, cough", use_llm?: boolean }
export async function predictDisease(req, res, next) {
  try {
    const userId = req.params.id || null;
    const { symptoms, use_llm = true } = req.body || {};

    if (!symptoms || !String(symptoms).trim()) {
      return res.status(400).json({ ok: false, error: 'symptoms is required' });
    }

    const result = await mcp.diagnose(symptoms, use_llm !== false);
    const data = extractToolJson(result);
    if (!data) {
      return res.status(502).json({ ok: false, error: 'Empty response from MCP diagnose tool', raw: result });
    }

    return res.status(200).json({ ok: true, userId, data });
  } catch (err) {
    next(err);
  }
}

// POST /predictTargetProtein
// Body: { symptoms: "fever, cough", max_targets?: number }
export async function predictTargetProtein(req, res, next) {
  try {
    const { symptoms, max_targets = 3 } = req.body || {};
    if (!symptoms || !String(symptoms).trim()) {
      return res.status(400).json({ ok: false, error: 'symptoms is required' });
    }

    const result = await mcp.analyzeTargets(symptoms, max_targets);
    const data = extractToolJson(result);
    if (!data) {
      return res.status(502).json({ ok: false, error: 'Empty response from MCP analyze_targets tool', raw: result });
    }

    // data shape: { target_protein_1: {...}, target_protein_2: {...}, ... }
    return res.status(200).json({ ok: true, targets: data });
  } catch (err) {
    next(err);
  }
}

// GET /getnewdrug/:id
// Query: ?name=EGFR&max=5&include_admet=true
export async function getnewdrug(req, res, next) {
  try {
    const protein_id = req.params.id || '';
    const {
      name: protein_name = '',
      max = 5,
      include_admet = 'true',
    } = req.query || {};

    if (!protein_id && !protein_name) {
      return res.status(400).json({ ok: false, error: 'Provide protein_id as path param or name as query (?name=)' });
    }

    const result = await mcp.ligandsForProtein({
      protein_id,
      protein_name,
      max_ligands: Number(max || 5),
      include_admet: String(include_admet).toLowerCase() !== 'false',
    });

    const data = extractToolJson(result);
    if (!data) {
      return res.status(502).json({ ok: false, error: 'Empty response from MCP get_ligands_for_protein tool', raw: result });
    }

    // data shape: { ligands: [...] }
    return res.status(200).json({ ok: true, protein_id, protein_name, ...data });
  } catch (err) {
    next(err);
  }
}

// GET /symptoms/:id
// Path param is a URL-encoded symptoms string. Example: /symptoms/headache%2C%20fever
export async function getSymptoms(req, res, next) {
  try {
    let encoded = req.params.id || '';
    let symptoms = '';
    try {
      symptoms = decodeURIComponent(encoded);
    } catch {
      symptoms = encoded;
    }

    const use_llm = String(req.query.use_llm || 'true').toLowerCase() !== 'false';

    if (!symptoms || !symptoms.trim()) {
      return res.status(400).json({ ok: false, error: 'Provide symptoms in the path as a URL-encoded string' });
    }

    const result = await mcp.diagnose(symptoms, use_llm);
    const data = extractToolJson(result);
    if (!data) {
      return res.status(502).json({ ok: false, error: 'Empty response from MCP diagnose tool', raw: result });
    }

    return res.status(200).json({ ok: true, input: symptoms, data });
  } catch (err) {
    next(err);
  }
}