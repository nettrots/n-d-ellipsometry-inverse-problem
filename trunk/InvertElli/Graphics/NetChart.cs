using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ChartDirector;


namespace Graphics
{
    public class NetChart:AbstractCharting
    {
        private double[] dataX;
        private double[] dataY;
        private double[] dataZ;
        private double min ;
        private double max ;

        private XYChart chart;
        private ContourLayer layer;
        private WinChartViewer viewer;

        public NetChart(WinChartViewer viewer)
        {
            this.viewer = viewer;
        }
        public double[] DataX
        {
            get { return dataX; }
            set { dataX = value; }
        }

        public double[] DataY
        {
            get { return dataY; }
            set { dataY = value; }
        }

        public double[] DataZ
        {
            get { return dataZ; }
            set { dataZ = value; }
        }
        public void ChangeData(List<double[]> data)
        {
            min = double.PositiveInfinity;
            max = double.NegativeInfinity;
            dataX=new double[data.Count];
            dataY = new double[data.Count];
            dataZ = new double[data.Count];
            for (int i = 0; i < data.Count; i++)
            {
                dataX[i] = data[i][0];
                dataY[i] = data[i][1];
                dataZ[i] = data[i][2];
                if (max < dataZ[i]) max = dataZ[i];
                if (min > dataZ[i]) min = dataZ[i];
            }
        }
        public override void Init()
        {
            chart = new XYChart(viewer.Width, viewer.Height, 0xffffff, 0x888888);
            chart.setSize(viewer.Width, viewer.Height);
            chart.setRoundedFrame();
            chart.setPlotArea(30, 30, viewer.Width - 60, viewer.Height - 60, 0);

        }

        public override void ChangeState()
        {
            Init();
            layer = chart.addContourLayer(dataX, dataY, dataZ);
        }

        public override void Draw()
        {
            
            chart.setClipping();
            chart.xAxis().setTickDensity(10);
            chart.yAxis().setTickDensity(10);
            // Move the grid lines in front of the contour layer
            chart.getPlotArea().moveGridBefore(layer);

            // Add a color axis (the legend) in which the top center is anchored at
            // (245, 455). Set the length to 330 pixels and the labels on the top
            // side.
            int w = chart.getWidth();
            int h = chart.getHeight();
            ColorAxis cAxis = layer.setColorAxis(0, 0, Chart.TopLeft, h - w - 40,
                Chart.Top);

            // Add a bounding box to the color axis using the default line color as
            // border.
            cAxis.setBoundingBox(Chart.Transparent, Chart.LineColor);


            // Set the color axis range as 0 to 20, with a step every 2 units
            cAxis.setLinearScale(min, max, Math.Abs(max - min) / 40);
            cAxis.setColorGradient();

            viewer.Chart= chart;
        }
    }
}
