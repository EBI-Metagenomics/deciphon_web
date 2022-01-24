import axios from "axios";
import config from "./config/config.json";

export const baseUrl = config.API_BASE;
export const pollingInterval = config.POLLING_INTERVAL_MS;

export default axios.create({
  baseURL: baseUrl,
});
