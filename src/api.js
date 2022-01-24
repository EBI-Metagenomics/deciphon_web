import axios from "axios";
import config from "./config/config.json";

export const baseUrl = config.API_BASE;

export default axios.create({
  baseURL: baseUrl,
});
