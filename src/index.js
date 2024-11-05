import React from "react";
import {
  createBrowserRouter,
  RouterProvider,
} from "react-router-dom";
import { createRoot } from 'react-dom/client';
import { ToastContainer } from "react-toastify";
import "react-toastify/dist/ReactToastify.css";
import "./App.css";
import Query from "./components/Query";
import Job, { loader } from "./components/Job";
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

const root = createRoot(document.getElementById("root"));

root.render(
  <React.StrictMode>
    <RouterProvider router={router} />
    <ToastContainer />
  </React.StrictMode>
);
