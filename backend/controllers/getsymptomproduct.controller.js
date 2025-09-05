import { PredictDisease } from "../models/newdrug.model.js";
import { ReactionResponse } from "../models/newdrug.model.js";
import { User } from "../models/auth.model.js";
import { MongoClient, ObjectId } from "mongodb";
import GeneratenewMolecule from "../models/generatenew.model.js";

export const getSymptomsAndAllProducts = async (req, res) => {
  try {
    const userId = req.params.id;

    // Verify user exists
    const user = await User.findById(userId);
    if (!user) {
      console.error(`User not found for userId: ${userId}`);
      return res.status(404).json({ success: false, message: "User not found" });
    }

    // MongoDB connection for reaction_responses
    const client = new MongoClient(process.env.MONGODB_URI);
    await client.connect();
    const db = client.db("test");
    const reactionResponsesCollection = db.collection("reaction_responses");

    // Fetch symptoms and productSmileId from PredictDisease
    const predictDiseases = await PredictDisease.find({ userId }).select("symptoms productSmileId");
    console.log("PredictDisease documents:", predictDiseases);

    if (!predictDiseases || predictDiseases.length === 0) {
      console.warn(`No symptoms found for userId: ${userId}`);
      await client.close();
      return res.status(404).json({ message: "No symptoms found for this user" });
    }

    // Form symptom groups dynamically from PredictDisease documents
    const symptomGroups = predictDiseases
      .map((doc) => [...new Set(doc.symptoms.filter(s => s && s.trim()))])
      .filter(group => group.length > 0);
    console.log("Symptom groups:", symptomGroups);

    if (symptomGroups.length === 0) {
      console.warn(`No valid symptom groups for userId: ${userId}`);
      await client.close();
      return res.status(400).json({ message: "No valid symptom groups found" });
    }

    // Collect productSmileIds
    const productSmileIds = predictDiseases
      .map((pd) => pd.productSmileId)
      .filter((id) => id)
      .map((id) => new ObjectId(id));
    console.log("Product Smile IDs:", productSmileIds);

    // Fetch product_smiles
    let productSmilesGroups = [];
    if (productSmileIds.length > 0) {
      const reactions = await reactionResponsesCollection
        .find({ _id: { $in: productSmileIds } }, { projection: { product_smiles: 1 } })
        .toArray();

      if (reactions.length > 0) {
        const allProductSmiles = reactions.map((r) => r.product_smiles.filter(s => s && s.trim()));
        if (symptomGroups.length === 1) {
          productSmilesGroups = [[...new Set(allProductSmiles.flat())]];
        } else {
          const groupSize = Math.ceil(allProductSmiles.length / symptomGroups.length);
          productSmilesGroups = symptomGroups.map((_, i) => {
            const start = i * groupSize;
            const end = start + groupSize;
            return [...new Set(allProductSmiles.slice(start, end).flat())];
          });
        }
      }
    }

    // Fallback to all product_smiles if needed
    if (productSmilesGroups.length < symptomGroups.length) {
      const allProducts = await reactionResponsesCollection
        .find({}, { projection: { product_smiles: 1 } })
        .toArray();
      const allProductSmiles = [
        ...new Set(allProducts.flatMap((doc) => doc.product_smiles).filter(s => s && s.trim())),
      ];

      if (allProductSmiles.length > 0) {
        const groupSize = Math.ceil(allProductSmiles.length / symptomGroups.length);
        productSmilesGroups = symptomGroups.map((_, i) => {
          const start = i * groupSize;
          const end = start + groupSize;
          return allProductSmiles.slice(start, end);
        });
      } else {
        productSmilesGroups = symptomGroups.map(() => []);
      }
    }
    console.log("Product SMILES groups before mapping:", productSmilesGroups);

    // Fetch molecules from GeneratenewMolecule to map SMILES to drug names
    const molecules = await GeneratenewMolecule.find({ userId }).select(
      "newSmiles originalSmiles newmoleculetitle potentialDiseases"
    );
    const smilesToDrugNameMap = {};
    const moleculeDetails = {};
    molecules.forEach((molecule) => {
      if (molecule.originalSmiles) {
        smilesToDrugNameMap[molecule.originalSmiles] = molecule.newSmiles;
      }
      if (!molecule.originalSmiles && molecule.newSmiles) {
        smilesToDrugNameMap[molecule.newSmiles] = molecule.newSmiles;
      }
      moleculeDetails[molecule.newSmiles] = {
        originalSmiles: molecule.originalSmiles,
        newmoleculetitle: molecule.newmoleculetitle,
        potentialDiseases: molecule.potentialDiseases,
      };
    });
    console.log("SMILES to drug name map:", smilesToDrugNameMap);

    // Map SMILES to drug names in productSmilesGroups
    productSmilesGroups = productSmilesGroups.map((group) =>
      group.map((smiles) => smilesToDrugNameMap[smiles] || smiles)
    );
    console.log("Product SMILES groups after mapping:", productSmilesGroups);

    // Close MongoDB connection
    await client.close();

    // Combine results
    const result = {
      symptoms: symptomGroups,
      productSmiles: productSmilesGroups,
      moleculeDetails,
      message: productSmilesGroups.every((group) => group.length === 0)
        ? "No product SMILES found in reaction responses"
        : undefined,
    };

    console.log(`Successfully fetched data for userId: ${userId}`, result);
    return res.status(200).json(result);
  } catch (error) {
    console.error("Error in getSymptomsAndAllProducts:", error);
    return res.status(500).json({
      message: "Server error",
      error: error.message,
    });
  }
};