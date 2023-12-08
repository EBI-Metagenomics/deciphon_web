import React from "react";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import Query from "./components/Query";
import Result from "./components/Result";
import About from "./components/About";

import { ToastContainer } from "react-toastify";
import "react-toastify/dist/ReactToastify.css";
import "./App.css";


export default function ParamsExample() {
  return (
    <Router>
      <div>
        <Routes>
          <Route path="/" element={<Query />} />
          <Route path="/about" element={<About />} />
          <Route path="/jobs/:jobid" element={<Result />} />
        </Routes>
        <ToastContainer />
      </div>
    </Router>
  );
}
