import mongoose from "mongoose";

export const connectionDb = async ()=>{
    try {
        const connect = await mongoose.connect(process.env.MONGODB_URI)
        console.log(`Mongodb connection: ${connect.connection.host}`)
        
    } catch (error) {
        console.log(`error in connection :${error.message}`)
        
    }

}