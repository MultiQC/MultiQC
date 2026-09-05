////////////////////////////////////////////////
// MultiQC Report Toolbox Utility Functions
////////////////////////////////////////////////

// Make functions available globally
window.hashCode = function (str) {
  var hash = 0;
  if (str.length == 0) return hash;
  for (i = 0; i < str.length; i++) {
    char = str.charCodeAt(i);
    hash = (hash << 5) - hash + char;
    hash = hash & hash; // Convert to 32bit integer
  }
  return hash;
};

// Convert RGB color to hex format
window.rgbToHex = function (rgb) {
  // Extract numbers from rgb(r, g, b) format
  const matches = rgb.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
  if (!matches) return rgb; // Return original if not RGB format

  // Convert each component to hex
  function componentToHex(c) {
    const hex = parseInt(c).toString(16);
    return hex.length === 1 ? "0" + hex : hex;
  }

  // Combine components with # prefix
  return "#" + componentToHex(matches[1]) + componentToHex(matches[2]) + componentToHex(matches[3]);
};

window.validate_regexp = function (pattern) {
  try {
    new RegExp(pattern, "g");
    return true;
  } catch (error) {
    showToast(
      "Invalid Regular Expression!",
      "Apologies, your regular expression pattern is invalid: <code>" +
        pattern +
        "</code><br><br>" +
        'For more help and testing, try it out at <a href="https://regex101.com/" target="_blank">regex101.com</a>.',
      "error",
    );
    return false;
  }
};

window.mqc_toolbox_confirmapply = function () {
  // Check if there's anything waiting to be applied
  if ($("#mqc_cols_apply").is(":enabled") && $("#mqc_cols").is(":visible")) {
    showToast(
      "Highlights Not Applied",
      "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      "warning",
    );
  }
  if ($("#mqc_renamesamples_apply").is(":enabled") && $("#mqc_renamesamples").is(":visible")) {
    showToast(
      "Rename Patterns Not Applied",
      "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      "warning",
    );
  }
  if ($("#mqc_hidesamples_apply").is(":enabled") && $("#mqc_hidesamples").is(":visible")) {
    showToast(
      "Hide Samples Not Applied",
      "Careful - your changes haven't been applied yet! Click the <em>Apply</em> button in the toolbox to set your changes.",
      "warning",
    );
  }
};

// Storage helper functions
window.saveToLocalStorage = function (key, value) {
  try {
    localStorage.setItem(key, value);
  } catch (e) {
    console.warn("Failed to save to localStorage:", e);
  }
};

window.getFromLocalStorage = function (key) {
  try {
    return localStorage.getItem(key);
  } catch (e) {
    console.warn("Failed to read from localStorage:", e);
    return null;
  }
};

window.notEmptyObj = function (obj) {
  return obj !== null && obj !== undefined && obj !== "";
};

window.dataUrlToBlob = function (dataUrl, mime) {
  // Split the data URL at the comma
  const byte_str = atob(dataUrl.split(",")[1]);
  const byte_numbers = new Array(byte_str.length);
  for (let i = 0; i < byte_str.length; i++) {
    byte_numbers[i] = byte_str.charCodeAt(i);
  }
  const byte_array = new Uint8Array(byte_numbers);
  return new Blob([byte_array], { type: mime });
};

window.addLogo = function (imageDataUrl, callback) {
  // Append "watermark" to the image
  let plotlyImage = new Image();
  plotlyImage.onload = function () {
    let canvas = document.createElement("canvas");
    let ctx = canvas.getContext("2d");

    // Set canvas size to double for retina display
    canvas.width = plotlyImage.width;
    canvas.height = plotlyImage.height; // additional height for the text

    ctx.drawImage(plotlyImage, 0, 0, plotlyImage.width, plotlyImage.height);

    // Text properties
    ctx.font = "9px system-ui"; // Set the desired font-size and type
    ctx.textAlign = "right";
    ctx.fillStyle = "#9f9f9f"; // Semi-transparent black text

    ctx.fillText("Created with MultiQC", plotlyImage.width - 15, plotlyImage.height - 6);
    // Callback with the combined image
    callback(canvas.toDataURL());
  };
  plotlyImage.src = imageDataUrl;
};
