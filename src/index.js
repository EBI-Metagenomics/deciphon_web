import React from "react";
import ReactDOM from "react-dom";
import {
  createBrowserRouter,
  RouterProvider,
} from "react-router-dom";
import { ToastContainer } from "react-toastify";
import "react-toastify/dist/ReactToastify.css";
import "./App.css";
import Query from "./components/Query";
import Job, {loader} from "./components/Job";
import About from "./components/About";

const router = createBrowserRouter([
  {
    path: "/",
    element: <Query />,
  },
  {
    path: "/about",
    element: <About />,
  },
  {
    path: "/jobs/:jobid",
    element: <Job />,
    loader: loader,
  },
]);

ReactDOM.render(
  <React.StrictMode>
    <RouterProvider router={router} />
    <ToastContainer />
  </React.StrictMode>,
  document.getElementById("root")
);
