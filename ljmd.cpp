#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#define NATOMS  1000 

using namespace std;

class Vector {
	//size definitions and length functions
	public:
		double x;
		double y;
		double z;
		Vector() {

		}
		Vector (double x_, double y_, double z_) {
			x = x_;
			y = y_;
			z = z_;
		}
		double len() {
			return sqrt(x*x + y*y + z*z);
		}
		double lenSqr() {
			return x*x + y*y + z*z;
		}
};

class Atom 
{
	//three vectors for pos/vel/force ->> specifies an atom
	public:
		Vector pos;
		Vector vel;
		Vector force;
		Atom (double x_, double y_, double z_) {
			pos = Vector(x_, y_, z_);
			vel = Vector(0, 0, 0);
			force = Vector(0, 0, 0);
		}
		// Kinetic energy for an atom given by .5*v^2
		double KE() {
			return .5 * vel.lenSqr();
		}
};

class BoxData 
// parameter definition for the box
{
	public:
		int natoms;				// number of atoms in box
		double lenWhole;		//full length	
		double lenHalf;			//half box length
		double vol;				//volume of the box
		double dens;			//density
		double temp;			//temperature of box
		double press;			//pressure in box
		double pressav;			//average pressure
		double cpress[6];		//critical pressure? but why 6 tems? -- related with the 6 elements in virial pressure def
		double pcor;			//pressure correction factor, arrise from the fact that atoms outside the cutoff radius contribute to pressure
		double delta_u;			//potential energy difference created by adding one particle to random position in box with no velocity
};

class SimData
//sinulation parameters defined
{
	public:
		int    ID;              //simulation id?
		double Tstar;			//T* = kBT/Є
		long   cyc_eq;			//number of equilibrium cycles
		long   cyc_pr;			//number of production steps/ cycles
		long   blockc;			//frequency for rescaling velocities
		long   blockd;			//frequency for blocking of data and output, i.e. time step between averaging and outputting data
		long   drift;			//frequency for checking momentum drift
		double rCut;			//cutoff radius	
		double rCutSqr;			//squared of cutoff	
		double dt;				//time step
};

class Energy 
//types of energies - variables defined
{
	public:
		double total;			//talk energy
		double poten;			//potential energy
		double kinet;			//kinetic energy
		double kinblock;		//?
		double potens;			//?
		double totals;			//total energy
};			

class Virial 
//virial coefficients for pressure equation, P = rho x kB x T* x (1 + B1 x rho + B2 x rho^2 + . . .)
{
	public:
		double total;		//
		double nbond[6];	// related with the boxdata cpress[6]
};

vector<Atom> atoms;
Virial pvir;
Energy en;
SimData sim;
BoxData box;

// declare functions that will be needed
void init_all();		
void read_sim();		
void write_config();	
void output();		
void enforcePBC();
void forces();
void integrate();
void kinet();
void initial_vel();
void initial_pos();
void outputXYZ(string fn);
void velscale();
void pressure();

// added code
void add_atom();
double poten_dif();

// declare random number generating functions
int  ran_num_int(int,int);
double  gauss(double);
double  ran_num_float(long ,int ,int);

int turn;

int main()
{
	FILE *sout; FILE *sout1; char name[50]; char name1[50];
	sprintf(name,"./simul.out");		
	sout = fopen(name,"w");
	fclose(sout);
	sprintf(name1,"./ener.out");
	sout1= fopen(name1,"w");
	fclose(sout1);
	ran_num_float(-5,0,1);	
	read_sim();		 
	init_all();		
	initial_pos();		  
	outputXYZ("./iniconfig.xyz");		 
	enforcePBC();			  
	initial_vel();		   
	kinet();			      
	forces();			      
	pressure();			  
	output();			      
	const int totalTurns = sim.cyc_eq + sim.cyc_pr;
	for (turn = 0; turn < totalTurns; turn++)
	{
		if (turn==sim.cyc_eq) {    
			fprintf(stdout,"Equilibration Completed\n");	
		}
		forces(); 		
		pressure();	

		// get potential energy here before the simulation integrates
		// we only want to run this part when the system has equilibrated
		if (turn>sim.cyc_eq) {
			double delta_u = poten_dif();
			box.delta_u = delta_u;	
		}	

		integrate();	
		kinet();			
		en.kinblock += en.kinet;		
		en.total = en.kinet + en.poten;
		enforcePBC();

		/* -------------------------------------------------------------- */
		/*                                                                */
		/* The following command does the isokinetic rescaling of         */
		/* velocities to achieve target temperature                       */
		/* Note in NVE simulation rescaling is done only during           */
		/* equilibration (i.e. when turn<sim.cyc_eq)                      */
		/*                                                                */
		/* -------------------------------------------------------------- */

		if (turn>1 && turn<sim.cyc_eq && (turn%sim.blockc ==0)) {
			velscale();
		}

		if ((turn% sim.blockd)==0 ) {
			output();
			outputXYZ("./config.xyz");
		}
	}
	return 0;
}

void read_sim ()
{
	char tt[80]; //for temp storage
	FILE *input;
	//read exiting files if present
	if ( NULL == (input=fopen("./simul.input","r")) ) {
		fprintf(stdout,"input file simul.input does not exist\n");
		exit(1);
	}
	input = fopen("./simul.input","r");
	sim.ID = 1;
	//fscanf reads inputs from the screen
	//format on fscanf  allows to read 1 significant digit after decimal
	//format on fgets -- pointer/length-to-read/stream-to-read
	fscanf(input, "%lf", &sim.Tstar);  fgets(tt, 80, input); 
	fscanf(input, "%ld", &sim.cyc_eq); fgets(tt, 80, input);
	fscanf(input, "%ld", &sim.cyc_pr); fgets(tt, 80, input);
	fscanf(input, "%ld", &sim.blockd); fgets(tt, 80, input);
	fscanf(input, "%ld", &sim.blockc); fgets(tt, 80, input);
	fscanf(input, "%ld", &sim.drift);  fgets(tt, 80, input);
	fscanf(input, "%lf", &sim.dt);     fgets(tt, 80, input);
	fscanf(input, "%lf", &sim.rCut);     fgets(tt, 80, input);
	fscanf(input, "%d",  &box.natoms); fgets(tt, 80, input);
	fscanf(input, "%lf", &box.dens);   fgets(tt, 80, input);
	fgets(tt,80,input);
	fclose(input);

}

void init_all ()
{
	// initialize all of the system variables
	box.temp = 0.0; 
	box.press = 0.0;
	box.pressav = 0.0;
	en.total = 0.0;
	en.poten = 0.0;
	en.kinet = 0.0;
	en.kinblock = 0.0;
	box.vol = (box.natoms)/(box.dens);
	box.lenWhole = pow(box.vol, 1.0/3); //full length = cuberoot of vol
	box.lenHalf = 0.5*box.lenWhole;
	sim.rCutSqr = sim.rCut * sim.rCut;
	box.pcor = (16.0/3.0)*M_PI*box.dens*box.dens*((2.0/3.0)*(pow(sim.rCut,-9.0))-(pow(sim.rCut, -3))); // pressure correction factor,
	// accounts for pressure resulting from particles outside the cutoff radius

	//box.delta_u = 0;
}

void initial_pos()
{
	int numPerDim = ceil(pow(box.natoms, 1.0/3));
	double boxLen = box.lenWhole;

	// loops over all occupied positions in a 3-D lattice
	// loop stops when the length of the vector atoms.size is the same as the preset box.natoms
	for (int xIdx=0; xIdx<numPerDim; xIdx++) {
		for (int yIdx=0; yIdx<numPerDim; yIdx++) {
			for (int zIdx=0; zIdx<numPerDim; zIdx++) {
				if ((int) atoms.size() == box.natoms) {
					return;
				} 
				// assigns x, y, and z values for atoms based on the iterative index and the box size
				double x = (xIdx + 0.5) * boxLen / numPerDim; //why add 0.5? -> places all of the atoms on one side of the box
				double y = (yIdx + 0.5) * boxLen / numPerDim;
				double z = (zIdx + 0.5) * boxLen / numPerDim;
				Atom a (x, y, z);
				atoms.push_back(a); // built in function to push another value onto vector set
			}
		}
	}
}

void initial_vel() {
	// assigns an initial value to each atom based on a Gaussian distribution of velocities,
	// then normalizes the velocities by subtracting the average, calculates an average kinetic energy in the box from that
	// and then normalizes the velocities based on the average kinetic energy and T*.
	// This ensures that the velocities are normally distributed, and are based on the input temperature.
	double sumvx = 0;
	double sumvy = 0;
	double sumvz = 0;
	double sumvSqr = 0;
	for (unsigned int i=0; i<atoms.size(); i++) {
		// iterates over total number of atoms, assigning each a velocity based on gaussian distribution
		Vector &vel = atoms[i].vel; // vel now references the velocity vector we created from each atom[i]
		vel.x = gauss(sim.Tstar);
		vel.y = gauss(sim.Tstar); 
		vel.z = gauss(sim.Tstar); 
		// running total velocity in each direction is calculated
		sumvx += vel.x;
		sumvy += vel.y;
		sumvz += vel.z;
	}
	for (unsigned int i=0; i<atoms.size(); i++) {
		// again iterating over all atoms in the box, average velocity in each direction is subtracted from the velocity of each atom
		Vector &vel = atoms[i].vel;
		vel.x -= sumvx/atoms.size(); // total velocity in x direction/ number of atoms = average x-dir velocity
		vel.y -= sumvy/atoms.size();
		vel.z -= sumvz/atoms.size();
		sumvSqr += 2 * atoms[i].KE(); // total v^2 = 2 x total kinetic energy, running total based on new normalized velocities
	}

	sumvSqr = sumvSqr/atoms.size(); // average v^2 = total v^2 / num of atoms
	for (unsigned int i=0; i<atoms.size(); i++) {
		// once more iterate over all atoms, and reassign velocities based on normalization (3/2 x T* / KE)^1/2
		Vector &vel = atoms[i].vel;
		vel.x *= sqrt(3*sim.Tstar/sumvSqr);
		vel.y *= sqrt(3*sim.Tstar/sumvSqr);
		vel.z *= sqrt(3*sim.Tstar/sumvSqr);
	}
}

void forces () {
	double enonbond  = 0.0;
	double enonbonds = 0.0;
	double wxx = 0.0;
	double wyx = 0.0;
	double wyy = 0.0;
	double wzx = 0.0;
	double wzy = 0.0;
	double wzz = 0.0;
	double boxLen, boxLenHalf, rCutSqr;
	double fx, fy, fz;
	double dx, dy, dz;
	double pbcx, pbcy, pbcz;
	double drSqr, sig=1, eps=1;
	double srSqr, sr6, sr12;
	double srSqrs, sr6s, sr12s, force;
	
	for (unsigned int i=0; i<atoms.size(); i++) {
		// reinitialize all forces to 0
		Vector &force = atoms[i].force; //using a reference here instead of a pointer.  Acts like a pointer, but syntax is like a normal object.
		force.x = 0;
		force.y = 0;
		force.z = 0;
	}

	boxLen = box.lenWhole;
	boxLenHalf = box.lenHalf;
	rCutSqr = sim.rCutSqr;
	const int numAtoms = atoms.size();
	for (int i=0; i<numAtoms-1; i++) {
		Atom &a = atoms[i];
		for (int j=i+1; j<numAtoms; j++) {
			Atom &b = atoms[j];

			dx = a.pos.x - b.pos.x;
			dy = a.pos.y - b.pos.y;
			dz = a.pos.z - b.pos.z;
			pbcx = 0.0;
			pbcy = 0.0;
			pbcz = 0.0;
			if(dx >  boxLenHalf) pbcx =- boxLen; // pbc -> periodic boundary condition
			if(dx < -boxLenHalf) pbcx =+ boxLen;
			if(dy >  boxLenHalf) pbcy =- boxLen;
			if(dy < -boxLenHalf) pbcy =+ boxLen;
			if(dz >  boxLenHalf) pbcz =- boxLen;
			if(dz < -boxLenHalf) pbcz =+ boxLen;
			dx += pbcx;
			dy += pbcy;
			dz += pbcz;
			drSqr = dx*dx + dy*dy + dz*dz;
			if (drSqr > rCutSqr) {
				continue;
			}
			srSqr  = (sig * sig) / drSqr;
			sr6  = srSqr * srSqr * srSqr;
			sr12 = sr6 * sr6;

			srSqrs  = (sig * sig) / rCutSqr;
			sr6s  = srSqrs * srSqrs * srSqrs;
			sr12s = sr6s * sr6s;

			force = (24.0*eps/drSqr) * (2.0*sr12 - sr6); // Lennard-Jones force
			enonbond  += 4.0 * eps * (sr12 - sr6);	// Lennard-Jones potential
			enonbonds += 4.0 * eps * (sr12s - sr6s);

			fx = force * dx;
			a.force.x += fx;
			b.force.x -= fx;
			fy = force * dy;
			a.force.y += fy;
			b.force.y -= fy;
			fz = force * dz;
			a.force.z += fz;
			b.force.z -= fz;

			wxx += fx * dx;
			wyx += fy * dx;
			wyy += fy * dy;
			wzx += fz * dx;
			wzy += fz * dy;
			wzz += fz * dz;
		}
	}
	// pvir is the first virial coefficient in P = rho x kB x T x (1 + B1xrho + B2xrho^2 + ....)
	pvir.nbond[0] = wxx;
	pvir.nbond[1] = wyx;
	pvir.nbond[2] = wyy;
	pvir.nbond[3] = wzx;
	pvir.nbond[4] = wzy;
	pvir.nbond[5] = wzz;
	pvir.total = (wxx+wyy+wzz)/3.0;

	en.poten = enonbond;
	en.potens = enonbond-enonbonds;
}

void integrate () {
	//get all the velocities from force results and positions
	double dt = sim.dt;

	for(unsigned int i=0; i<atoms.size(); i++) {
		Atom &a = atoms[i];
		a.vel.x += dt * a.force.x;
		a.vel.y += dt * a.force.y;
		a.vel.z += dt * a.force.z;
		
		a.pos.x += dt * a.vel.x;
		a.pos.y += dt * a.vel.y;
		a.pos.z += dt * a.vel.z;

	}
}

void output () {
	//store in file "simul.out":  iterations, temp, av_pres, tmep
	//? should we not be storing the time evolution of the system?
	//what is turn?
	FILE *sout;FILE *sout1; char name[50]; char name1[50];
	sprintf(name,"./simul.out");
	sout = fopen(name,"a");
	box.pressav /= sim.blockd; 
	fprintf(stdout,"%d iterations completed, pressure = %lf, temperature = %lf \n",turn, box.pressav, box.temp);
	fprintf(sout,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",turn,box.temp,box.press,en.poten,en.kinet,en.total,box.dens,box.delta_u);
	//fprintf(sout," Cycle / Time step:             %12ld    %12.4f \n",turn,tt);
	//fprintf(sout," Temperature:                   %12.2f \n",box.temp);
	//fprintf(sout," Density:                       %12.4f \n",box.dens);
	//fprintf(sout," Length:                        %12.4f \n",box.lenWhole);
	//fprintf(sout," Volume:                        %12.4f \n",box.vol);
	//fprintf(sout," Pressure (atomic/molecular):   %12.4f \n",box.press);
	//fprintf(sout," Press.(xx,yy,zz):              %12.4f  %12.4f  %12.4f \n",box.cpress[0],box.cpress[2],box.cpress[5]);
	//fprintf(sout," Press.(xy,xz,yz):              %12.4f  %12.4f  %12.4f \n",box.cpress[1],box.cpress[3],box.cpress[4]);
	//fprintf(stdout," Potential Energy:              %12.4f \n",en.poten);
	//fprintf(sout," Kinetic Energy:                %12.4f \n",en.kinet);
	//fprintf(sout," Total Energy:                  %12.4f \n",en.total);

	fprintf(stdout, "potential difference = %lf \n", box.delta_u); // test pot dif

	sprintf(name1,"./ener.out");
	sout1= fopen(name1,"a");
	fprintf(sout1,"%d\t%lf\t%lf\t%lf\n",turn,en.kinet,en.poten,en.total); //store KE PE E_T in a different file
	fclose (sout1);
	fclose(sout);
	box.pressav = 0.0; // resets average pressure
}

void kinet () {//get total KE and temperature from there
	en.kinet = 0.0;
	for (unsigned int i=0; i<atoms.size(); i++) {
		Atom &a = atoms[i];
		en.kinet += a.KE(); // en.kinet is the total kinetic energy of all atoms
	}
	box.temp = 2.0*en.kinet/(3.0*atoms.size()); // T = (2/3)*KE_average, KE_average = KE_total/(# of atoms)
}

void enforcePBC () {
	//this function moves the atoms that have exited the left wall back in from the right wall and so on
	double boxLen = box.lenWhole;
	double boxLenHalf = box.lenHalf;
	double pbcx,pbcy,pbcz;

	for (unsigned int i=0; i<atoms.size(); i++) {
		Vector &pos = atoms[i].pos;

		pbcx = 0.0;
		pbcy = 0.0;
		pbcz = 0.0;

		if(pos.x >  boxLenHalf) pbcx =- boxLen;
		if(pos.x < -boxLenHalf) pbcx =+ boxLen;
		if(pos.y >  boxLenHalf) pbcy =- boxLen;
		if(pos.y < -boxLenHalf) pbcy =+ boxLen;
		if(pos.z >  boxLenHalf) pbcz =- boxLen;
		if(pos.z < -boxLenHalf) pbcz =+ boxLen;

		pos.x += pbcx;
		pos.y += pbcy;
		pos.z += pbcz;

	}
}

void outputXYZ (string fn) {
	char t[2] ="C";
	FILE *io;
	//io=fopen("./iniconfig.xyz","w");
	io=fopen(fn.c_str(), "w");
	fprintf(io,"   %d\n",box.natoms);
	fprintf(io,"\n");
	for (unsigned int i=0; i<atoms.size(); i++) {
		fprintf(io," %2s   ", t);
		fprintf(io," %15.6lf  ", atoms[i].pos.x);
		fprintf(io," %15.6lf  ", atoms[i].pos.y);
		fprintf(io," %15.6lf\n", atoms[i].pos.z);
	}
	fclose(io);
}

void velscale() {
	// scales the velocities back to a gaussian distribution based on T* and kinetic energies
	en.kinblock /= (sim.blockc*box.natoms); // averages kinetic energy
	for (unsigned int i=0; i<atoms.size(); i++) {
		Atom &a = atoms[i];
		a.vel.x *= sqrt(3*sim.Tstar/(2*en.kinblock));
		a.vel.y *= sqrt(3*sim.Tstar/(2*en.kinblock));
		a.vel.z *= sqrt(3*sim.Tstar/(2*en.kinblock));
	}
	en.kinblock = 0.0; // resets kinetic energy counter to zero
}

void pressure() {
	// calculates pressure
	box.press = box.dens * box.temp + pvir.total/box.vol + box.pcor; // P* = rho* x T* + B1/V* + correction factor
	box.pressav += (box.press);
	box.cpress[0] = box.dens * box.temp + pvir.nbond[0]/box.vol+ box.pcor;
	box.cpress[1] = box.dens * box.temp + pvir.nbond[1]/box.vol+ box.pcor;
	box.cpress[2] = box.dens * box.temp + pvir.nbond[2]/box.vol+ box.pcor;
	box.cpress[3] = box.dens * box.temp + pvir.nbond[3]/box.vol+ box.pcor;
	box.cpress[4] = box.dens * box.temp + pvir.nbond[4]/box.vol+ box.pcor;
	box.cpress[5] = box.dens * box.temp + pvir.nbond[5]/box.vol+ box.pcor;
}

// new function to add atom to the box with ranodom position and no velocity
void add_atom() {

	int numPerDim = ceil(pow(box.natoms, 1.0/3)); // defines number of lattice sites occupied, not sure about this line
	double boxLen = box.lenWhole;
	// Get a random position for the new atom.  We will assume that the probability of placing it directly on an existing atom is very small,
	// and for computational times' sake will not check that this is true.
	double xran = ran_num_int(0,numPerDim);
	double yran = ran_num_int(0,numPerDim);
	double zran = ran_num_int(0,numPerDim);
	double xpos = (xran + 0.5) * boxLen / numPerDim;
	double ypos = (yran + 0.5) * boxLen / numPerDim;
	double zpos = (zran + 0.5) * boxLen / numPerDim;

	Atom a(xpos, ypos, zpos);
	atoms.push_back(a);
};

double poten_dif() {

	// obtains a value for potential difference, u1 - u0 by adding an atom at a random position, calculating potential energy in the box,
	// and subtracting the potential energy when the atom is removed again.
	add_atom(); // add the random atom
	forces(); // calculate forces
	double u1 = en.poten; // get potential energy with extra atom
	atoms.pop_back(); // built in function to remove the last item in a vector, here the randomly placed atom
	forces(); // recalculate forces
	double u0 = en.poten; // get potential energy without extra atom
	double delta_u = u1 - u0; // take the difference in potential energy
	//fprintf(stdout,"%f,%f,%f \n",u1,u0,delta_u);

	return delta_u;
};



/* ##### Random number generator code ####### */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 0.0
#define RNMX 1.0

//random number on the gaussian distro from ones on the normal distr. 
//implemented here is the  polar form of the Box-Muller transformation
double gauss(double sigma) {
	double ran1,ran2,ransq,v;
	do {
		ran1=2.0*ran_num_float(1,0,1)-1.0; 
		ran2=2.0*ran_num_float(1,0,1)-1.0;
		ransq=ran1*ran1+ran2*ran2;
	} while(ransq >= 1.0);

	v=sqrt(sigma)*ran1*sqrt(-2.0*log(ransq)/ransq);
	return v;
}

// Function to generate floating point random numbers.

double ran_num_float(long idum,int range1,int range2) {
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (idum <= 0) {
		if (-(idum) < 1) {
			idum=1;
		} else {
			idum = -(idum);
		}

		idum2=(idum);

		for (j=NTAB+7;j>=0;j--)	{
			k=(idum)/IQ1;
			idum=IA1*(idum-k*IQ1)-IR1*k;
			if (idum < 0) idum += IM1;
			if (j < NTAB) iv[j] = idum;
		}
		iy=iv[0];
	}

	k=(idum)/IQ1;
	idum=IA1*(idum-k*IQ1)-IR1*k;

	if (idum < 0) {
		idum += IM1;
	}
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;

	if (idum2 < 0) {
		idum2+=IM2;
	}

	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j]=idum;
	if (iy<1) {
		iy+=IMM1;
	}
	temp=(double)AM*iy;

	return range1+(range2-range1)*temp;
}

// Function to generate integer type random numbers.

int ran_num_int(int range1,int range2) {
	return (int)ran_num_float(1,range1,range2);
}

/* #### END OF RANDOM NUMBER GENERATOR CODE ########## */
