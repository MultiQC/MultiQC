////////////////////////////////////////////////
// Call all deferred plotting functions
////////////////////////////////////////////////

$(function () {
  decompressPlotData(mqc_compressed_plotdata, (mqc_plotdata, err) => {
    if (err) {
      console.error(err);
      return;
    }
    callAfterDecompressed.forEach(function (fn) {
      fn(mqc_plotdata);
    });
  });
});
