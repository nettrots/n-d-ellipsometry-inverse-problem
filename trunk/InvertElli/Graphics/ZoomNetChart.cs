using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using ChartDirector;
using System.Drawing;

namespace Graphics
{
    public class ZoomNetChart:AbstractCharting
    {
        private XYChart chart;
        private WinChartViewer viewer;
        private Image old_Image;
        private Image new_Image;
        private System.Drawing.Graphics g;
        public ZoomNetChart(WinChartViewer viewer, AbstractCharting basecharter)
        {
            this.viewer = viewer;
            old_Image = viewer.Image;
            g=System.Drawing.Graphics.FromImage(new_Image=(Image)old_Image.Clone());

            
        }
        public override void Init()
        {
            
        }

        public enum ClickState
        {
            FIRST_POINT,LAST_POINT,CLEAR_POINTS
        }

        private int X1, Y1, X2, Y2;

        public override void ChangeState(params object[] args)
        {
            switch ((ClickState)args[0])
            {
                case ClickState.FIRST_POINT:
                    X1 = ((Point) args[1]).X;
                    Y1 = ((Point) args[1]).Y;
                    break;
                case ClickState.LAST_POINT:
                    X2 = ((Point) args[1]).X;
                    Y2 = ((Point) args[1]).Y;
                    break;
                case  ClickState.CLEAR_POINTS:
                    X1 = 0;
                    X2 = 0;
                    Y1 = 0;
                    Y1 = 0;
                    //originalDrawArea.
                    break;
            }
        }

        public override void Draw()
        {
            if (!(Y2 == 0 && X2 == 0))
            {
                g.DrawRectangle(new Pen(Color.Black,2),new Rectangle(X1,Y1,X1+X2,Y1+Y2) );
            }
            else
            {
                new_Image = (Image)old_Image.Clone();
            }

            viewer.Image = new_Image;
        }
    }
}
