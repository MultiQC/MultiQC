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
        // TODO: Should be able to remove this as soon as we use Bootstrap toasts instead of jquery.toast.css
        assetFileNames: (assetInfo) => {
          if (assetInfo.name.endsWith(".css")) {
            return "css/multiqc.min.css";
          }
          return "assets/[name][extname]";
        },
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
