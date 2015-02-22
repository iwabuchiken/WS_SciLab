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