/*!
 * Color mode toggler for MultiQC based on Bootstrap's implementation
 * Licensed under the Apache License, Version 2.0
 */

(() => {
  "use strict";

  // Check if dark mode is enabled via config
  if (window.mqc_config && window.mqc_config.template_dark_mode === false) {
    return;
  }

  const getStoredTheme = () => localStorage.getItem("mqc-theme");
  const setStoredTheme = (theme) => localStorage.setItem("mqc-theme", theme);

  const getPreferredTheme = () => {
    const storedTheme = getStoredTheme();
    if (storedTheme && storedTheme !== "auto") {
      return storedTheme;
    }

    // Only use system preference for "auto" or no stored theme
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  };

  const setTheme = (theme) => {
    if (theme === "auto") {
      document.documentElement.setAttribute(
        "data-bs-theme",
        window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light",
      );
    } else {
      document.documentElement.setAttribute("data-bs-theme", theme);
    }
  };

  setTheme(getPreferredTheme());

  const showActiveTheme = (theme, focus = false) => {
    const themeSwitcher = document.querySelector("#bd-theme");

    if (!themeSwitcher) {
      return;
    }

    const themeSwitcherText = document.querySelector("#bd-theme-text");
    const activeThemeIcon = document.querySelector(".theme-icon-active");
    const btnToActive = document.querySelector(`[data-bs-theme-value="${theme}"]`);

    if (!btnToActive) {
      return;
    }

    // Update all dropdown items
    document.querySelectorAll("[data-bs-theme-value]").forEach((element) => {
      element.classList.remove("active");
      element.setAttribute("aria-pressed", "false");
      // Hide all checkmarks
      const checkIcon = element.querySelector(".check-icon");
      if (checkIcon) {
        checkIcon.classList.add("d-none");
      }
    });

    // Activate the selected theme
    btnToActive.classList.add("active");
    btnToActive.setAttribute("aria-pressed", "true");

    // Show checkmark for active theme
    const activeCheckIcon = btnToActive.querySelector(".check-icon");
    if (activeCheckIcon) {
      activeCheckIcon.classList.remove("d-none");
    }

    // Update the main toggle button icon by copying from the active button
    if (activeThemeIcon && btnToActive) {
      const activeIcon = btnToActive.querySelector(".me-2");
      if (activeIcon) {
        activeThemeIcon.innerHTML = activeIcon.innerHTML;
      }
    }

    const themeSwitcherLabel = `${themeSwitcherText.textContent} (${btnToActive.dataset.bsThemeValue})`;
    themeSwitcher.setAttribute("aria-label", themeSwitcherLabel);

    if (focus) {
      themeSwitcher.focus();
    }
  };

  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => {
    const storedTheme = getStoredTheme();
    // Only update theme if user has "auto" selected or no theme stored
    if (storedTheme === "auto" || !storedTheme) {
      setTheme(getPreferredTheme());
    }
  });

  window.addEventListener("DOMContentLoaded", () => {
    showActiveTheme(getPreferredTheme());

    document.querySelectorAll("[data-bs-theme-value]").forEach((toggle) => {
      toggle.addEventListener("click", () => {
        const theme = toggle.getAttribute("data-bs-theme-value");
        setStoredTheme(theme);
        setTheme(theme);
        showActiveTheme(theme, true);
      });
    });
  });
})();
