module.exports = {
  printWidth: 120,
  tabWidth: 2,
  trailingComma: "all",
  plugins: [require.resolve("prettier-plugin-jinja-template")],
  overrides: [
    {
      files: ["*.html"],
      options: {
        parser: "jinja-template",
      },
    },
  ],
};
