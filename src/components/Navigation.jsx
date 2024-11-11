export default function Navigation({ page }) {
  return (
    <nav className="vf-navigation vf-navigation--on-this-page vf-u-fullbleed vf-cluster">
      <ul className="vf-navigation__list vf-list vf-cluster__inner">
        <li className="vf-navigation__item">
          <a href="/" className="vf-navigation__link" aria-current={page === "query" && "page"}>Query</a>
        </li>
        <li className="vf-navigation__item">
          <a href="/about" className="vf-navigation__link" aria-current={page === "about" && "page"}>About</a>
        </li>
      </ul>
    </nav>
  );
}
