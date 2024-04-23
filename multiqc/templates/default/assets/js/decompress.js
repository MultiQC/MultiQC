function base64ToBytes(base64) {
  const binString = atob(base64);
  return Uint8Array.from(binString, (m) => m.codePointAt(0));
}

async function decompress_mqc_plotdata(data) {
  const data_bytes = base64ToBytes(data);
  const data_blob = new Blob([data_bytes]);
  const decompressor = new DecompressionStream("gzip");
  const decompressed_stream = data_blob.stream().pipeThrough(decompressor);
  const decompressed_bytes = await new Response(decompressed_stream).blob();
  const encoded_text_buffer = await decompressed_bytes.arrayBuffer();
  const text_decoder = new TextDecoder();
  text_decoder.encoding = "utf-8";
  const decoded_text = text_decoder.decode(encoded_text_buffer);
  return decoded_text;
}
