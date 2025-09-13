# src/rag/main.py
import os
import json
from pathlib import Path
from typing import List, Dict, Tuple, Any, Optional

import numpy as np
import pandas as pd
from openai import OpenAI

# ---------- Config via env ----------
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY", "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")
OPENROUTER_BASE_URL = os.getenv("OPENROUTER_BASE_URL", "https://openrouter.ai/api/v1")
CHAT_MODEL = os.getenv("OPENROUTER_CHAT_MODEL", "deepseek/deepseek-chat-v3.1:free")
# Important: use vendor-prefixed name on OpenRouter
EMBEDDING_MODEL = os.getenv("OPENROUTER_EMBED_MODEL", "openai/text-embedding-3-small")

# ---------- Helpers ----------
def normalize_symptom(s: str) -> str:
    if s is None:
        return ""
    s = s.strip().lower()
    s = s.replace("-", " ").replace("_", " ")
    s = " ".join(s.split())
    return s

def normalize_symptom_list(symptoms: List[str]) -> Tuple[str, ...]:
    cleaned = [normalize_symptom(s) for s in symptoms if s and str(s).strip()]
    cleaned = sorted(list({c for c in cleaned if c}))
    return tuple(cleaned)

def parse_symptom_text(symptom_text: str) -> List[str]:
    if not symptom_text:
        return []
    parts = [p.strip() for p in symptom_text.replace(";", ",").split(",")]
    return [p for p in parts if p]

def default_dataset_paths() -> Tuple[str, str]:
    base = Path(__file__).parent / "Datasets"
    return str(base / "DiseaseAndSymptoms.csv"), str(base / "DiseasePrecaution.csv")

# ---------- Data loading ----------
def load_datasets(symptoms_csv: str, precautions_csv: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    return pd.read_csv(symptoms_csv), pd.read_csv(precautions_csv)

def build_symptom_index(ds: pd.DataFrame) -> Dict[Tuple[str, ...], List[str]]:
    disease_col = next((c for c in ds.columns if c.lower() in ("disease", "diseases", "label")), None)
    if not disease_col:
        raise ValueError("No 'Disease' column found in DiseaseAndSymptoms.csv")

    symptom_cols = [c for c in ds.columns if c != disease_col and ("symptom" in c.lower() or c.lower() == "symptoms")]
    index: Dict[Tuple[str, ...], List[str]] = {}

    for _, row in ds.iterrows():
        disease = str(row[disease_col]).strip()
        raw_symptoms: List[str] = []
        # If a single 'Symptoms' column exists
        symptoms_col_name = next((c for c in ds.columns if c.lower() == "symptoms"), None)
        if symptoms_col_name:
            cell = row.get(symptoms_col_name, "")
            if pd.isna(cell):
                cell = ""
            raw_symptoms = [s.strip() for s in str(cell).split(",") if s and str(s).strip()]
        else:
            # Multiple symptom columns
            for c in symptom_cols:
                val = row.get(c, "")
                if pd.isna(val) or str(val).strip() == "":
                    continue
                raw_symptoms.append(str(val).strip())

        combo = normalize_symptom_list(raw_symptoms)
        if not combo:
            continue
        index.setdefault(combo, []).append(disease)

    return index

def build_precautions_map(dp: pd.DataFrame) -> Dict[str, List[str]]:
    disease_col = next((c for c in dp.columns if c.lower() in ("disease", "diseases", "label")), None)
    if not disease_col:
        raise ValueError("No 'Disease' column found in DiseasePrecaution.csv")

    precaution_cols = [c for c in dp.columns if "precaution" in c.lower()]
    if not precaution_cols:
        raise ValueError("No precaution columns found in DiseasePrecaution.csv")

    prec_map: Dict[str, List[str]] = {}
    for _, row in dp.iterrows():
        disease = str(row[disease_col]).strip()
        vals: List[str] = []
        for c in precaution_cols:
            v = row.get(c, "")
            if pd.isna(v) or str(v).strip() in ("", "nan"):
                continue
            vals.append(str(v).strip())
        prec_map[disease] = vals
    return prec_map

# ---------- OpenRouter ----------
def make_openrouter_client() -> OpenAI:
    if not OPENROUTER_API_KEY:
        raise RuntimeError("OPENROUTER_API_KEY not set")
    return OpenAI(base_url=OPENROUTER_BASE_URL, api_key=OPENROUTER_API_KEY)

def embed_texts(client: OpenAI, texts: List[str], model: Optional[str] = None) -> np.ndarray:
    model = model or EMBEDDING_MODEL

    # Call embeddings endpoint directly
    resp = client.embeddings.create(model=model, input=texts)

    # --- Fix: handle string or dict response from OpenRouter ---
    if isinstance(resp, str):  # OpenRouter sometimes returns JSON as string
        try:
            resp = json.loads(resp)
        except json.JSONDecodeError:
            raise RuntimeError(f"Unexpected string response from embeddings: {resp}")

    if isinstance(resp, dict):  # OpenRouter native JSON
        data = resp.get("data", [])
    else:  # OpenAI SDK object
        data = getattr(resp, "data", [])

    # Convert embeddings into numpy
    vecs = [np.array(item["embedding"], dtype=np.float32) for item in data]
    return np.vstack(vecs) if vecs else np.zeros((0, 0), dtype=np.float32)


def chat_complete(client: OpenAI, system_prompt: str, user_prompt: str, model: Optional[str] = None) -> str:
    resp = client.chat.completions.create(
        model=model or CHAT_MODEL,
        messages=[{"role": "system", "content": system_prompt}, {"role": "user", "content": user_prompt}],
        temperature=0.2,
    )
    # --- Fix: handle if resp is str or dict ---
    if isinstance(resp, str):
        try:
            resp = json.loads(resp)
        except json.JSONDecodeError:
            raise RuntimeError(f"Unexpected string response from chat: {resp}")
    
    if isinstance(resp, dict):
        return resp.get("choices", [{}])[0].get("message", {}).get("content", "")
    else:
        return resp.choices[0].message.content or ""

# ---------- Simple vector store ----------
class SimpleVectorStore:
    def __init__(self, dim: int):
        self.dim = dim
        self._vecs = np.empty((0, dim), dtype=np.float32)
        self._meta: List[Dict[str, Any]] = []

    @staticmethod
    def _normalize(mat: np.ndarray) -> np.ndarray:
        norms = np.linalg.norm(mat, axis=1, keepdims=True) + 1e-12
        return mat / norms

    def add(self, vectors: np.ndarray, metadatas: List[Dict[str, Any]]):
        if vectors.shape[1] != self.dim or vectors.shape[0] != len(metadatas):
            raise ValueError("Vector dims/length mismatch")
        v = self._normalize(vectors.astype(np.float32))
        self._vecs = np.vstack([self._vecs, v]) if self._vecs.size else v
        self._meta.extend(metadatas)

    def search(self, query_vec: np.ndarray, top_k: int = 5):
        q = self._normalize(query_vec.astype(np.float32).reshape(1, -1))
        sims = (self._vecs @ q.T).ravel()
        idxs = np.argsort(-sims)[:top_k]
        return [(float(sims[i]), self._meta[i]) for i in idxs]

def build_documents(index: Dict[Tuple[str, ...], List[str]]):
    texts, metas = [], []
    for combo, diseases in index.items():
        for d in diseases:
            texts.append(f"Disease: {d}\nSymptoms: {'; '.join(combo)}")
            metas.append({"disease": d, "symptoms": list(combo)})
    return texts, metas

def prepare_vector_store(client: OpenAI, index: Dict[Tuple[str, ...], List[str]]) -> Optional[SimpleVectorStore]:
    texts, metas = build_documents(index)
    if not texts:
        return None
    vecs = embed_texts(client, texts)
    if vecs.size == 0:
        return None
    store = SimpleVectorStore(dim=vecs.shape[1])
    store.add(vecs, metas)
    return store

# ---------- Core logic ----------
def find_exact_match(user_symptoms: List[str], symptom_index: Dict[Tuple[str, ...], List[str]]):
    combo = normalize_symptom_list(user_symptoms)
    diseases = symptom_index.get(combo)
    if diseases:
        return diseases[0], combo
    return None, combo

def get_precautions(disease: str, prec_map: Dict[str, List[str]]) -> List[str]:
    return prec_map.get(disease, [])

def generate_llm_answer(client: OpenAI, disease: str, precautions: List[str], user_symptoms: List[str]) -> str:
    sys = "Parse the vector store that contain that csv so here I want if the combination of that symptoms exists for only single disease then only give me 100% match otherwise no and don't predict means like don't give outout like flu etc., noting"
    usr = (
        f"User symptoms: {', '.join(user_symptoms)}\n"
        f"Matched disease: {disease}\n"
        "Precautions:\n" + "\n".join(f"- {p}" for p in precautions) + "\n"
        "Write a short, clear response summarizing the disease and the precautions."
    )
    try:
        return chat_complete(client, sys, usr, CHAT_MODEL)
    except Exception:
        if not precautions:
            return f"Disease: {disease}\nNo explicit precautions found."
        return "Disease: {}\nPrecautions:\n{}".format(disease, "\n".join(f"- {p}" for p in precautions))

# ---------- Orchestration ----------
def run_pipeline(
    user_symptoms_input: str,
    symptoms_csv: Optional[str] = None,
    precautions_csv: Optional[str] = None,
    use_llm: bool = True,
) -> Dict[str, Any]:
    if not user_symptoms_input or not str(user_symptoms_input).strip():
        return {"status": "error", "message": "No symptoms provided."}

    if symptoms_csv is None or precautions_csv is None:
        def_sym, def_prec = default_dataset_paths()
        symptoms_csv = symptoms_csv or def_sym
        precautions_csv = precautions_csv or def_prec

    user_symptoms = parse_symptom_text(user_symptoms_input)

    ds, dp = load_datasets(symptoms_csv, precautions_csv)
    index = build_symptom_index(ds)
    prec_map = build_precautions_map(dp)

    client = make_openrouter_client()
    # _ = prepare_vector_store(client, index)  # built but not used for matching per your spec
    # Removed the above line to avoid error, as OpenRouter does not support embeddings API and this is unused.

    matched_disease, normalized_combo = find_exact_match(user_symptoms, index)
    if matched_disease:
        precautions = get_precautions(matched_disease, prec_map)
        answer = generate_llm_answer(client, matched_disease, precautions, list(normalized_combo)) if use_llm \
                 else ("Disease: {}\n{}".format(matched_disease, "\n".join(f"- {p}" for p in precautions)))
        return {
            "status": "ok",
            "match_type": "exact",
            "normalized_symptoms": list(normalized_combo),
            "disease": matched_disease,
            "precautions": precautions,
            "answer": answer,
        }

    return {
        "status": "new_combination",
        "normalized_symptoms": list(normalized_combo),
        "message": "symptoms combination is new it is new dieasese",
    }

