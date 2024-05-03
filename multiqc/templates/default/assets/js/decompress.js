function decompressPlotData(base64Str, callback) {
  // Decode the Base64 string to bytes.
  const binaryString = atob(base64Str);
  const bytes = Uint8Array.from(binaryString, (m) => m.codePointAt(0));

  // Check if DecompressionStream is supported
  if ("DecompressionStream" in window) {
    // Passing a callback to work around that DecompressionStream is async
    decompressUsingDecompressionStream(bytes, callback);
  } else {
    decompressUsingPako(bytes, callback);
  }
}

function decodeDecompressedBytes(decompressedBytes) {
  const decoder = new TextDecoder("utf-8");
  const jsonStr = decoder.decode(decompressedBytes);

  let data = JSON.parse(jsonStr);
  return data;
}

function decompressUsingDecompressionStream(bytes, callback) {
  const ds = new DecompressionStream("gzip");
  const decompressedStream = new Response(bytes.buffer).body.pipeThrough(ds);

  new Response(decompressedStream)
    .arrayBuffer()
    .then((decompressedArrayBuffer) => {
      const decompressedBytes = new Uint8Array(decompressedArrayBuffer);
      callback(decodeDecompressedBytes(decompressedBytes), null);
    })
    .catch((error) => {
      console.error("Decompression with DecompressionStream failed:", error);
      return decompressUsingPako(bytes, callback);
    });
}

function decompressUsingPako(bytes, callback) {
  try {
    const decompressedBytes = pako.inflate(bytes);
    callback(decodeDecompressedBytes(decompressedBytes), null);
  } catch (error) {
    console.error("Decompression with pako failed:", error);
    callback(null, error); // Error callback
  }
}
