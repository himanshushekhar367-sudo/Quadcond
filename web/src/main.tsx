import { StrictMode } from "react";
import { createRoot } from "react-dom/client";
import { AppShell } from "@/components/layout/AppShell";
import "@/styles.css";

const root = document.getElementById("root");
if (!root) throw new Error("no #root element in index.html");

createRoot(root).render(
  <StrictMode>
    <AppShell />
  </StrictMode>,
);
