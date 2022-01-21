import { render, screen } from "@testing-library/react";
import App from "./App";

test("renders query page", () => {
  render(<App />);
  const titleElement = screen.getByText(/query deciphon/i);
  expect(titleElement).toBeInTheDocument();
});
