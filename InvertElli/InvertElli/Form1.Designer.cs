namespace InvertElli
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.button1 = new System.Windows.Forms.Button();
            this.textBox1 = new System.Windows.Forms.TextBox();
            this.delta_tb = new System.Windows.Forms.TextBox();
            this.psi_tb = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.dstep_b = new System.Windows.Forms.TextBox();
            this.nstep_b = new System.Windows.Forms.TextBox();
            this.nmin_b = new System.Windows.Forms.TextBox();
            this.nmax_b = new System.Windows.Forms.TextBox();
            this.dmin_b = new System.Windows.Forms.TextBox();
            this.dmax_b = new System.Windows.Forms.TextBox();
            this.label9 = new System.Windows.Forms.Label();
            this.label10 = new System.Windows.Forms.Label();
            this.textBox2 = new System.Windows.Forms.TextBox();
            this.amb_n = new System.Windows.Forms.TextBox();
            this.label11 = new System.Windows.Forms.Label();
            this.l2_n = new System.Windows.Forms.TextBox();
            this.l1_n = new System.Windows.Forms.TextBox();
            this.subst_n = new System.Windows.Forms.TextBox();
            this.textBox7 = new System.Windows.Forms.TextBox();
            this.l1_d = new System.Windows.Forms.TextBox();
            this.l2_d = new System.Windows.Forms.TextBox();
            this.textBox10 = new System.Windows.Forms.TextBox();
            this.textBox11 = new System.Windows.Forms.TextBox();
            this.label12 = new System.Windows.Forms.Label();
            this.label13 = new System.Windows.Forms.Label();
            this.label14 = new System.Windows.Forms.Label();
            this.label15 = new System.Windows.Forms.Label();
            this.label16 = new System.Windows.Forms.Label();
            this.label17 = new System.Windows.Forms.Label();
            this.label18 = new System.Windows.Forms.Label();
            this.done_l = new System.Windows.Forms.Label();
            this.label20 = new System.Windows.Forms.Label();
            this.subst_k = new System.Windows.Forms.TextBox();
            this.l1_k = new System.Windows.Forms.TextBox();
            this.l2_k = new System.Windows.Forms.TextBox();
            this.textBox15 = new System.Windows.Forms.TextBox();
            this.amb_k = new System.Windows.Forms.TextBox();
            this.label19 = new System.Windows.Forms.Label();
            this.aoi = new System.Windows.Forms.TextBox();
            this.SuspendLayout();
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(475, 30);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(75, 23);
            this.button1.TabIndex = 0;
            this.button1.Text = "run";
            this.button1.UseVisualStyleBackColor = true;
            this.button1.Click += new System.EventHandler(this.button1_Click);
            // 
            // textBox1
            // 
            this.textBox1.Anchor = ((System.Windows.Forms.AnchorStyles)(((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom)
                        | System.Windows.Forms.AnchorStyles.Left)));
            this.textBox1.Location = new System.Drawing.Point(1, 12);
            this.textBox1.Multiline = true;
            this.textBox1.Name = "textBox1";
            this.textBox1.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBox1.Size = new System.Drawing.Size(410, 336);
            this.textBox1.TabIndex = 1;
            // 
            // delta_tb
            // 
            this.delta_tb.Location = new System.Drawing.Point(450, 126);
            this.delta_tb.Name = "delta_tb";
            this.delta_tb.Size = new System.Drawing.Size(100, 20);
            this.delta_tb.TabIndex = 2;
            this.delta_tb.Text = "160";
            // 
            // psi_tb
            // 
            this.psi_tb.Location = new System.Drawing.Point(450, 153);
            this.psi_tb.Name = "psi_tb";
            this.psi_tb.Size = new System.Drawing.Size(100, 20);
            this.psi_tb.TabIndex = 3;
            this.psi_tb.Text = "17";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(557, 126);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(30, 13);
            this.label1.TabIndex = 5;
            this.label1.Text = "delta";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(557, 160);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(20, 13);
            this.label2.TabIndex = 6;
            this.label2.Text = "psi";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(417, 199);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(29, 13);
            this.label3.TabIndex = 9;
            this.label3.Text = "nmin";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(417, 225);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(32, 13);
            this.label4.TabIndex = 10;
            this.label4.Text = "nmax";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(417, 308);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(32, 13);
            this.label5.TabIndex = 14;
            this.label5.Text = "dmax";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(417, 283);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(29, 13);
            this.label6.TabIndex = 13;
            this.label6.Text = "dmin";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(569, 308);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(39, 13);
            this.label7.TabIndex = 18;
            this.label7.Text = "d_step";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(569, 218);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(39, 13);
            this.label8.TabIndex = 17;
            this.label8.Text = "n_step";
            // 
            // dstep_b
            // 
            this.dstep_b.Location = new System.Drawing.Point(519, 301);
            this.dstep_b.Name = "dstep_b";
            this.dstep_b.Size = new System.Drawing.Size(44, 20);
            this.dstep_b.TabIndex = 16;
            this.dstep_b.Text = "0.5";
            // 
            // nstep_b
            // 
            this.nstep_b.Location = new System.Drawing.Point(519, 218);
            this.nstep_b.Name = "nstep_b";
            this.nstep_b.Size = new System.Drawing.Size(44, 20);
            this.nstep_b.TabIndex = 15;
            this.nstep_b.Text = "0.05";
            // 
            // nmin_b
            // 
            this.nmin_b.Location = new System.Drawing.Point(452, 199);
            this.nmin_b.Name = "nmin_b";
            this.nmin_b.Size = new System.Drawing.Size(44, 20);
            this.nmin_b.TabIndex = 19;
            this.nmin_b.Text = "1.4";
            // 
            // nmax_b
            // 
            this.nmax_b.Location = new System.Drawing.Point(452, 222);
            this.nmax_b.Name = "nmax_b";
            this.nmax_b.Size = new System.Drawing.Size(44, 20);
            this.nmax_b.TabIndex = 20;
            this.nmax_b.Text = "1.6";
            // 
            // dmin_b
            // 
            this.dmin_b.Location = new System.Drawing.Point(450, 280);
            this.dmin_b.Name = "dmin_b";
            this.dmin_b.Size = new System.Drawing.Size(44, 20);
            this.dmin_b.TabIndex = 21;
            this.dmin_b.Text = "50";
            // 
            // dmax_b
            // 
            this.dmax_b.Location = new System.Drawing.Point(450, 305);
            this.dmax_b.Name = "dmax_b";
            this.dmax_b.Size = new System.Drawing.Size(44, 20);
            this.dmax_b.TabIndex = 22;
            this.dmax_b.Text = "200";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(417, 264);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(70, 13);
            this.label9.TabIndex = 23;
            this.label9.Text = "thicknessin A";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(417, 183);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(70, 13);
            this.label10.TabIndex = 24;
            this.label10.Text = "index number";
            // 
            // textBox2
            // 
            this.textBox2.Location = new System.Drawing.Point(677, 224);
            this.textBox2.Name = "textBox2";
            this.textBox2.ReadOnly = true;
            this.textBox2.Size = new System.Drawing.Size(44, 20);
            this.textBox2.TabIndex = 26;
            this.textBox2.Text = "x";
            // 
            // amb_n
            // 
            this.amb_n.Location = new System.Drawing.Point(677, 199);
            this.amb_n.Name = "amb_n";
            this.amb_n.Size = new System.Drawing.Size(44, 20);
            this.amb_n.TabIndex = 25;
            this.amb_n.Text = "1.333";
            this.amb_n.TextChanged += new System.EventHandler(this.anb_n_TextChanged);
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(697, 158);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(60, 13);
            this.label11.TabIndex = 27;
            this.label11.Text = "layer model";
            // 
            // l2_n
            // 
            this.l2_n.Location = new System.Drawing.Point(677, 250);
            this.l2_n.Name = "l2_n";
            this.l2_n.Size = new System.Drawing.Size(44, 20);
            this.l2_n.TabIndex = 28;
            this.l2_n.Text = "1.5";
            // 
            // l1_n
            // 
            this.l1_n.Location = new System.Drawing.Point(677, 276);
            this.l1_n.Name = "l1_n";
            this.l1_n.Size = new System.Drawing.Size(44, 20);
            this.l1_n.TabIndex = 29;
            this.l1_n.Text = "1.4599";
            // 
            // subst_n
            // 
            this.subst_n.Location = new System.Drawing.Point(677, 303);
            this.subst_n.Name = "subst_n";
            this.subst_n.Size = new System.Drawing.Size(44, 20);
            this.subst_n.TabIndex = 30;
            this.subst_n.Text = "3.8858";
            // 
            // textBox7
            // 
            this.textBox7.Location = new System.Drawing.Point(773, 303);
            this.textBox7.Name = "textBox7";
            this.textBox7.ReadOnly = true;
            this.textBox7.Size = new System.Drawing.Size(44, 20);
            this.textBox7.TabIndex = 35;
            this.textBox7.Text = "x";
            // 
            // l1_d
            // 
            this.l1_d.Location = new System.Drawing.Point(773, 276);
            this.l1_d.Name = "l1_d";
            this.l1_d.Size = new System.Drawing.Size(44, 20);
            this.l1_d.TabIndex = 34;
            this.l1_d.Text = "10";
            // 
            // l2_d
            // 
            this.l2_d.Location = new System.Drawing.Point(773, 250);
            this.l2_d.Name = "l2_d";
            this.l2_d.Size = new System.Drawing.Size(44, 20);
            this.l2_d.TabIndex = 33;
            this.l2_d.Text = "30";
            // 
            // textBox10
            // 
            this.textBox10.Location = new System.Drawing.Point(773, 224);
            this.textBox10.Name = "textBox10";
            this.textBox10.ReadOnly = true;
            this.textBox10.Size = new System.Drawing.Size(44, 20);
            this.textBox10.TabIndex = 32;
            this.textBox10.Text = "x";
            // 
            // textBox11
            // 
            this.textBox11.Location = new System.Drawing.Point(773, 199);
            this.textBox11.Name = "textBox11";
            this.textBox11.ReadOnly = true;
            this.textBox11.Size = new System.Drawing.Size(44, 20);
            this.textBox11.TabIndex = 31;
            this.textBox11.Text = "x";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(632, 199);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(44, 13);
            this.label12.TabIndex = 36;
            this.label12.Text = "ambient";
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(632, 224);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(25, 13);
            this.label13.TabIndex = 37;
            this.label13.Text = "Film";
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(632, 250);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(38, 13);
            this.label14.TabIndex = 39;
            this.label14.Text = "layer 2";
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(632, 276);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(38, 13);
            this.label15.TabIndex = 38;
            this.label15.Text = "layer 1";
            // 
            // label16
            // 
            this.label16.AutoSize = true;
            this.label16.Location = new System.Drawing.Point(674, 176);
            this.label16.Name = "label16";
            this.label16.Size = new System.Drawing.Size(13, 13);
            this.label16.TabIndex = 41;
            this.label16.Text = "n";
            // 
            // label17
            // 
            this.label17.AutoSize = true;
            this.label17.Location = new System.Drawing.Point(626, 306);
            this.label17.Name = "label17";
            this.label17.Size = new System.Drawing.Size(50, 13);
            this.label17.TabIndex = 40;
            this.label17.Text = "substract";
            // 
            // label18
            // 
            this.label18.AutoSize = true;
            this.label18.Location = new System.Drawing.Point(770, 176);
            this.label18.Name = "label18";
            this.label18.Size = new System.Drawing.Size(34, 13);
            this.label18.TabIndex = 43;
            this.label18.Text = "d in A";
            // 
            // done_l
            // 
            this.done_l.AutoSize = true;
            this.done_l.Location = new System.Drawing.Point(684, 30);
            this.done_l.Name = "done_l";
            this.done_l.Size = new System.Drawing.Size(0, 13);
            this.done_l.TabIndex = 44;
            // 
            // label20
            // 
            this.label20.AutoSize = true;
            this.label20.Location = new System.Drawing.Point(720, 176);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(13, 13);
            this.label20.TabIndex = 50;
            this.label20.Text = "k";
            // 
            // subst_k
            // 
            this.subst_k.Location = new System.Drawing.Point(723, 303);
            this.subst_k.Name = "subst_k";
            this.subst_k.Size = new System.Drawing.Size(44, 20);
            this.subst_k.TabIndex = 49;
            this.subst_k.Text = "-0.02";
            // 
            // l1_k
            // 
            this.l1_k.Location = new System.Drawing.Point(723, 276);
            this.l1_k.Name = "l1_k";
            this.l1_k.Size = new System.Drawing.Size(44, 20);
            this.l1_k.TabIndex = 48;
            this.l1_k.Text = "0";
            this.l1_k.TextChanged += new System.EventHandler(this.textBox13_TextChanged);
            // 
            // l2_k
            // 
            this.l2_k.Location = new System.Drawing.Point(723, 250);
            this.l2_k.Name = "l2_k";
            this.l2_k.Size = new System.Drawing.Size(44, 20);
            this.l2_k.TabIndex = 47;
            this.l2_k.Text = "0";
            // 
            // textBox15
            // 
            this.textBox15.Location = new System.Drawing.Point(723, 224);
            this.textBox15.Name = "textBox15";
            this.textBox15.ReadOnly = true;
            this.textBox15.Size = new System.Drawing.Size(44, 20);
            this.textBox15.TabIndex = 46;
            this.textBox15.Text = "x";
            // 
            // amb_k
            // 
            this.amb_k.Location = new System.Drawing.Point(723, 199);
            this.amb_k.Name = "amb_k";
            this.amb_k.Size = new System.Drawing.Size(44, 20);
            this.amb_k.TabIndex = 45;
            this.amb_k.Text = "0";
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(643, 129);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(33, 13);
            this.label19.TabIndex = 52;
            this.label19.Text = "angle";
            // 
            // aoi
            // 
            this.aoi.Location = new System.Drawing.Point(687, 126);
            this.aoi.Name = "aoi";
            this.aoi.Size = new System.Drawing.Size(44, 20);
            this.aoi.TabIndex = 51;
            this.aoi.Text = "60";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(829, 378);
            this.Controls.Add(this.label19);
            this.Controls.Add(this.aoi);
            this.Controls.Add(this.label20);
            this.Controls.Add(this.subst_k);
            this.Controls.Add(this.l1_k);
            this.Controls.Add(this.l2_k);
            this.Controls.Add(this.textBox15);
            this.Controls.Add(this.amb_k);
            this.Controls.Add(this.done_l);
            this.Controls.Add(this.label18);
            this.Controls.Add(this.label16);
            this.Controls.Add(this.label17);
            this.Controls.Add(this.label14);
            this.Controls.Add(this.label15);
            this.Controls.Add(this.label13);
            this.Controls.Add(this.label12);
            this.Controls.Add(this.textBox7);
            this.Controls.Add(this.l1_d);
            this.Controls.Add(this.l2_d);
            this.Controls.Add(this.textBox10);
            this.Controls.Add(this.textBox11);
            this.Controls.Add(this.subst_n);
            this.Controls.Add(this.l1_n);
            this.Controls.Add(this.l2_n);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.textBox2);
            this.Controls.Add(this.amb_n);
            this.Controls.Add(this.label10);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.dmax_b);
            this.Controls.Add(this.dmin_b);
            this.Controls.Add(this.nmax_b);
            this.Controls.Add(this.nmin_b);
            this.Controls.Add(this.label7);
            this.Controls.Add(this.label8);
            this.Controls.Add(this.dstep_b);
            this.Controls.Add(this.nstep_b);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.psi_tb);
            this.Controls.Add(this.delta_tb);
            this.Controls.Add(this.textBox1);
            this.Controls.Add(this.button1);
            this.Name = "Form1";
            this.Text = "Ellipsometry";
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.TextBox textBox1;
        private System.Windows.Forms.TextBox delta_tb;
        private System.Windows.Forms.TextBox psi_tb;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.TextBox dstep_b;
        private System.Windows.Forms.TextBox nstep_b;
        private System.Windows.Forms.TextBox nmin_b;
        private System.Windows.Forms.TextBox nmax_b;
        private System.Windows.Forms.TextBox dmin_b;
        private System.Windows.Forms.TextBox dmax_b;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.TextBox textBox2;
        private System.Windows.Forms.TextBox amb_n;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.TextBox l2_n;
        private System.Windows.Forms.TextBox l1_n;
        private System.Windows.Forms.TextBox subst_n;
        private System.Windows.Forms.TextBox textBox7;
        private System.Windows.Forms.TextBox l1_d;
        private System.Windows.Forms.TextBox l2_d;
        private System.Windows.Forms.TextBox textBox10;
        private System.Windows.Forms.TextBox textBox11;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.Label label16;
        private System.Windows.Forms.Label label17;
        private System.Windows.Forms.Label label18;
        private System.Windows.Forms.Label done_l;
        private System.Windows.Forms.Label label20;
        private System.Windows.Forms.TextBox subst_k;
        private System.Windows.Forms.TextBox l1_k;
        private System.Windows.Forms.TextBox l2_k;
        private System.Windows.Forms.TextBox textBox15;
        private System.Windows.Forms.TextBox amb_k;
        private System.Windows.Forms.Label label19;
        private System.Windows.Forms.TextBox aoi;
    }
}

