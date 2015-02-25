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