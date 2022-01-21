import React from "react";
import { BrowserRouter as Router, Routes, Route } from "react-router-dom";
import Query from "./components/Query";
import Result from "./components/Result";
import { ToastContainer } from "react-toastify";
import "react-toastify/dist/ReactToastify.css";
import "./App.css";

export default function ParamsExample() {
  return (
    <Router>
      <div>
        <Routes>
          <Route path="/" element={<Query />} />
          <Route path="/results/:jobid" element={<Result />} />
        </Routes>
        <ToastContainer />
      </div>
    </Router>
  );
}
