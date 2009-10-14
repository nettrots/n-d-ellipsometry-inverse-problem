using System;
using System.Collections.Generic;
using System.Drawing;
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
        private void ChangeData(List<double[]> data)
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

        private int marginX=50;
        private int marginY = 30;
        private int marginXright = 4 * 30;
        private int marginYbottom = 2 * 30;
        public override void Init()
        {
            chart = new XYChart(viewer.Width, viewer.Height, 0xffffff, 0x888888);
            chart.setSize(viewer.Width, viewer.Height);
            chart.setRoundedFrame();
            chart.setPlotArea(marginX, marginY, viewer.Width - marginXright, viewer.Height -  marginYbottom, Chart.metalColor(Color.Gray));
            chart.xAxis().setTickDensity(10);
            chart.yAxis().setTickDensity(10);
         
        }

        public PointF getXYcoord(int x,int y)
        {
            x = x - marginX;
            y = y - marginY;
            if (x < 0 || x > viewer.Width - marginXright + marginX) return new PointF();
            if (y < 0 || y > viewer.Height - marginYbottom + marginY) return new PointF();
            float newX = (x) / (float)(viewer.Width - marginXright );
            float newY = 1 - (y) / (float)(viewer.Height - marginYbottom );
            xmn = (float)chart.xAxis().getMinValue();
            xmx = (float)chart.xAxis().getMaxValue();
            ymn = (float)chart.yAxis().getMinValue();
            ymx = (float)chart.yAxis().getMaxValue();
           
            return new PointF((xmx - xmn) * newX + xmn, (ymx - ymn) * newY + ymn);
        }
        float xmn;
        float ymn;
        float xmx;
        float ymx;

        public override void ChangeState(params object[] arguments)
        {
            Init();
            ChangeData((List<double[]>) arguments[0]);
            layer = chart.addContourLayer(dataX, dataY, dataZ);
            
        }

        public override void Draw()
        {
            
            // Move the grid lines in front of the contour layer
            //chart.getPlotArea().moveGridBefore(layer);

            // Add a color axis (the legend) in which the top center is anchored at
            // (245, 455). Set the length to 330 pixels and the labels on the top
            // side.
            int w = chart.getWidth();
            int h = chart.getHeight();

//            ColorAxis cAxis = layer.setColorAxis(w-3*margin, 0, Chart.TopLeft2, w/2,
//                Chart.Left);
//            cAxis.setBoxMargin(2);
//            cAxis.setBoundingBox(Chart.Transparent);
//            // Set the color axis range as 0 to 20, with a step every 2 units
//            //cAxis.setTickDensity(1,-1);
//            cAxis.setLinearScale(min, max, Math.Abs(max - min) / 10);
//            cAxis.setColorGradient();

            chart.layout();
            chart.layoutAxes();
            viewer.Chart = chart;
           
        }
    }
}
