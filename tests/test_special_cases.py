# def test_special_cases(data_dir, tmp_path):
#     multiqc.run(
#         testing.data_dir() / "special_cases",
#         cfg=ClConfig(
#             strict=True,
#             output_dir=tmp_path,
#         ),
#     )
#     assert len(report.general_stats_data) > 0 or sum(len(m.sections) for m in report.modules) > 0
