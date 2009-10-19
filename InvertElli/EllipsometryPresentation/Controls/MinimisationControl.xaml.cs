using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using ChartDirector;

namespace EllipsometryPresentation.Controls
{
    /// <summary>
    /// Interaction logic for MinimisationMethodComponent.xaml
    /// </summary>
    public partial class MinimisationControl: UserControl
    {
        private WinChartViewer viewer;
        public MinimisationControl()
        {
            // Create WinForms Usercontrol and 
            // add it to the WindowsFormsHost
            
            InitializeComponent();
            viewer = new WinChartViewer();
            windowsFormsHost.Child = viewer;
            viewer.BackgroundImage =
                Bitmap.FromFile("C:/Documents and Settings/Admin/Мои документы/Мои рисунки/Ux1.bmp");

        }
    }
}
