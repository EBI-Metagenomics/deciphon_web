import axios from "axios";
import config from "./config/config.json";

export default axios.create({
  baseURL: config.API_BASE,
});
