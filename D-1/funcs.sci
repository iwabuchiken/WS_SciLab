//	plotcomplex
//	funcgrid
//	plot_3d
//	plot_3d_sphere
//	plot_3d_sphere_v2
//	plot_3d_Exp_Sin
//	plot_Exp_Sin
//	mesh_grid
abc = 'abc'

function []=plotcomplex()

	//REF how to code http://www.openeering.com/sites/default/files/Plotting_in_Scilab.pdf
	//REF how to plot http://www.matrixlab-examples.com/complex-numbers.html
	t = linspace(0,1,101);
	A = exp(%i*t); 
	y = imag(A);
	x = real(A);
	
	plot(x,y);

endfunction

function []=funcgrid()

	s=poly(0,'s'); //Defines s to be a polynomial variable
	K=1; T=1; //Gain and time-constant
	sys=syslin('c',K/(T*s+1)); //Creates sys as a continuous-time ('c') syslin model. 
	t=[0:0.05:5]; //Time vector to be used in simulation
	y1=csim('step',t,sys); //Simulates system sys with step input
	scf(1);clf; //Opens and clears figure 1
	plot(t,y1)
	ax1=gca();ax1.grid=[0,0]; //Adds grid to the plot
	u=sin(5*t); //Creates input signal
	y2=csim(u,t,sys); //Simulates system sys with u as input
	scf(2);clf;
	plot(t,y2);
	ax1=gca();ax1.grid=[0,0];

endfunction

function []=plot_3d()
	
	t=[0:0.3:2*%pi]';
	z=sin(t)*cos(t');
	plot3d(t,t,z)	

endfunction

function []=plot_3d_sphere()
	
	x = -1:0.1:1;
	//y = 1-sqrt(x);
	y = sqrt(1-x.^2);
	
	z = sqrt(1-x.^2-y.^2);
	
	plot3d(x,y,z);
	
	//z=sin(t)*cos(t');
	//plot3d(t,t,z)	

endfunction

function []=plot_3d_sphere_v2()
	
	a = -0.9:0.1:0.9;
	b = a;
	
	c = sqrt(1-a.^2-b.^2);
	//z = sqrt(1-x.^2-y.^2);
	
	plot3d(a,b,c);
	//plot3d(x,y,z);
	
	//z=sin(t)*cos(t');
	//plot3d(t,t,z)	

endfunction

function []=plot_3d_Exp_Sin()

	x = -10:0.1:10;
	y = 0:0.1:10;
	
	[xc yc] = meshgrid(x, y);
	w = %pi/4;
	b = -0.5;
	
	plot(x,sin(x));

	scf(2);
	
	x = -5:0.1:5;
	y = 0:0.1:10;
	[xc yc] = meshgrid(x, y);
	w = %pi/4;
	b = -0.5;
	z = exp(b*yc).*sin(w*xc);
	plot3d(x, y, z);
		
endfunction

function [x, y, xc, yc]=plot_Exp_Sin()

	x = -1:0.1:1;
	y = 0:0.1:1;
	
	[xc yc] = meshgrid(x, y);

	return [x, y, xc, yc];		
	
endfunction

function [x,y,xc,yc,yc_e,yc_e_t,yc_e_t1, yc_t, yc_t1]=mesh_grid()

	x = -1:0.1:1;
	y = 0:0.1:1;
	
	[xc yc] = meshgrid(x, y);

	// exponentials
	yc_e = exp(yc);
	
	yc_e_t = yc_e';
	yc_e_t1 = yc_e_t(1:1,1:11);	//=> range of: row 1 to 1, column 1 to 11
	
	// transpose
	yc_t = yc';
	yc_t1 = yc_t(1:1,1:11);
	
	// plot
	plot(yc_t1, yc_e_t1);
	
	// return
	return [x,y,xc,yc,yc_e,yc_e_t,yc_e_t1, yc_t, yc_t1];		
	
endfunction

function [x, y]=palabora()

	//x = -1:0.1:1;
	x = 10:0.01:11;

	y = sqrt(x - 10);
	//y = sqrt(x - 10);
	
	return [x, y];		
	//return [x, y, xc, yc];		
	
endfunction

function [x1, y1, x2,y2]=palabora2()

	x1_left = 3;
	x1_right = 4;

	//x = -1:0.1:1;
	range1 = x1_left:0.01:x1_right;
	range2 = -11:0.01:-10;

	y1 = sqrt(range1 - x1_left);
	y2 = -sqrt(range1 - x1_left);
	
	return [range1, y1 range1, y2];		
	
endfunction

//function [x1, y1, x2,y2]=palabora3()
function [rng1, y1, rng1, y2, rng2, y3, rng2, y4]=palabora3()

	x1_left = 3;
	x1_right = 4;

	x2_left = -4;
	x2_right = -3;

	rng1 = x1_left:0.01:x1_right;

	rng2 = x2_left:0.01:x2_right;

	// 1st, 2nd
	y1 = sqrt(rng1.^2 - x1_left^2);
	y2 = -sqrt(rng1.^2 - x1_left^2);
	
	// 3rd, 4th
	y3 = sqrt(rng2.^2 - x2_right^2);
	y4 = -sqrt(rng2.^2 - x2_right^2);
	
	return [rng1, y1, rng1, y2, rng2, y3, rng2, y4];		
	
endfunction

function [rng1, y1, rng1, y2, rng2, y3, rng2, y4]=palabora4()

	x1_left = 3;
	x1_right = 4;

	x2_left = -4;
	x2_right = -3;

	rng1 = x1_left:0.01:x1_right;

	rng2 = x2_left:0.01:x2_right;

	// 1st, 2nd
	y1 = (1/2) * sqrt((1/2*rng1).^2 - (1/2*x1_left)^2);
	y2 = -(1/2)*sqrt((1/2*rng1).^2 - (1/2*x1_left)^2);
	//y1 = (1/2) * sqrt(rng1.^2 - x1_left^2);
	//y2 = -(1/2)*sqrt(rng1.^2 - x1_left^2);
	
	// 3rd, 4th
	y3 = sqrt(rng2.^2 - x2_right^2);
	y4 = -sqrt(rng2.^2 - x2_right^2);
	
	return [rng1, y1, rng1, y2, rng2, y3, rng2, y4];		
	
endfunction

function [x1, x2, rng1]=palabora_crossY()

	t1_top = 4;
	t1_bottom = 3;

	rng1 = t1_top:-0.01:t1_bottom;

	// 1st, 2nd
	x1 = sqrt((1/2*rng1).^2 - (1/2*t1_bottom)^2);
	x2 = -sqrt((1/2*rng1).^2 - (1/2*t1_bottom)^2);
	
	return [x1, x2, rng1];		
	
endfunction

// @return => (t1-t0)/(t'1-t'0)
function [b, a]=lorenz_trans_time()

	b = 0:0.1:0.9;
	
	//a = 1 / sqrt(1-b.^2);
	a = 1 ./ sqrt(1-b.^2);

	return [b, a];		
	//return [b, a'];		
	
endfunction

function []=func()

	//return [x, y, xc, yc];		
	
endfunction

