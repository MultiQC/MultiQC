import { resolve } from "path";

export default {
  root: resolve(__dirname, "src"),
  build: {
    outDir: "../compiled",
    rollupOptions: {
      input: {
        main: resolve(__dirname, "src/js/main.js"),
      },
      output: {
        entryFileNames: "js/multiqc.min.js",
        assetFileNames: "css/multiqc.min.css",
      },
    },
    minify: "terser",
    cssMinify: true,
  },
  css: {
    preprocessorOptions: {
      scss: {
        silenceDeprecations: ["import", "mixed-decls", "color-functions", "global-builtin"],
      },
    },
  },
};
