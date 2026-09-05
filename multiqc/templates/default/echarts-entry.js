// Spec build + VisualMapComponent + CanvasRenderer.
// VisualMap is not in the plan's component list but the prototype heatmap
// (02, correlation matrix) uses it for its colour legend, so measure it too.
import * as echarts from "echarts/core";
import { BarChart, LineChart, ScatterChart, CustomChart, HeatmapChart, BoxplotChart } from "echarts/charts";
import {
  TooltipComponent,
  LegendComponent,
  GridComponent,
  DataZoomComponent,
  MarkLineComponent,
  MarkAreaComponent,
  VisualMapComponent,
  TitleComponent,
  ToolboxComponent,
} from "echarts/components";
import { SVGRenderer, CanvasRenderer } from "echarts/renderers";

echarts.use([
  BarChart,
  LineChart,
  ScatterChart,
  CustomChart,
  HeatmapChart,
  BoxplotChart,
  TooltipComponent,
  LegendComponent,
  GridComponent,
  DataZoomComponent,
  MarkLineComponent,
  MarkAreaComponent,
  VisualMapComponent,
  TitleComponent,
  // Plotly-style click+drag box zoom (POLISH.md #17): the toolbox's `dataZoom` feature
  // is what powers ECharts' rubber-band box-zoom cursor mode (`takeGlobalCursor`, see
  // echarts-plotting.js), even though the toolbox's own icon row is kept off-canvas.
  ToolboxComponent,
  SVGRenderer,
  CanvasRenderer,
]);

export * from "echarts/core";
export { echarts };
