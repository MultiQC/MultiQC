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
  SVGRenderer,
  CanvasRenderer,
]);

export * from "echarts/core";
export { echarts };
