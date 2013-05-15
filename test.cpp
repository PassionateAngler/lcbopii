#include <gmtl/gmtl.h>
#include <gmtl/Matrix.h>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <deque>
#include <cmath>
#include "atom.h"
#include "lcbopii.h"

using namespace  simul;

void test_G(LCBOPII lcbopii)
{
	   double y = -1.0;
	   double step =2/50.0; // 1 - (-1)
	   double z = 0.0;
	   double y0;

	   while(z <= 8)
	   {
		   y = -1.0;
		   std::cout << "# z=" << z << std::endl;

		   while(y <= 1.0)
		   {
			   std::cout << y << " " << lcbopii.G(y,z) << std::endl;
			   //std::cout << y << " " << lcbopii.G_1(y) << std::endl;
			   //std::cout << std::pow((1 - y), 2) << std::endl;
			   y += step;
		   }
		   y0 = lcbopii.y_0(z);
		   //std::cout << "G_1(y0) = " << lcbopii.G_1(y0) << " G_2(y0) = " << lcbopii.G_2(y0, z, y0) << std::endl;
		   //std::cout << "G_1_prim(y0) = " << lcbopii.G_1_prim(y0)
			//	     << " G_2_prim(y0) = " << lcbopii.G_2_prim(y0, z, y0) << std::endl;
		   std::cout << std::endl;
		   std::cout << std::endl;
		   z += 1.0;
	   }
}

void test_H(LCBOPII lcbopii)
{
	double x = -1.25;
	double step = (2.5)/400;

	while(x < 1.25)
	{

		//std::cout << "# H_2(d) =" << lcbopii.H_2(lcbopii.d) << std::endl;
		//std::cout << "# H_2_prim(d) =" << lcbopii.H_2_prim(lcbopii.d) << std::endl;
		//std::cout << "# H_3(d) =" << lcbopii.H_3(lcbopii.d) << std::endl;
		std::cout << x << " " << lcbopii.H(x) << std::endl;
		x += step;
	}
}

int main()
{
   //Atom::bond_type atoms, atoms2;
   /*std::deque<Atom *> atoms;
   for(int i=0; i<10; i++)
   {

	   Atom * a = new Atom(i*2.0, i/2.0, i/3.0);
	   a->id = i;
	   atoms.push_back(a);
   }

   for(std::deque<Atom *>::iterator it = atoms.begin(); it != atoms.end()-1; it++)
   {
	   Atom * a = *it;
	   a->addBond(*(it+1));
   }

   atoms[3]->removeBond(atoms[0]);

   for(std::deque<Atom *>::iterator it = atoms.begin(); it != atoms.end(); it++)
   {
	   Atom * a = *it;
	   std::cout<< a->id << " bonds:" << std::endl;
	   for(Atom::bond_type::iterator jt = a->get_bonds().begin(); jt != a->get_bonds().end(); jt++)
	   {
		   std::cout << "\t" << (*jt)->id << std::endl;
	   }
   }*/

   LCBOPII lcbopii;// = new LCBOPII();

   test_H(lcbopii);
   /*std::cout << lcbopii.H(lcbopii.d) << std::endl;
   std::cout << "H_2(d)" << lcbopii.H_2(lcbopii.d) << std::endl;
   std::cout << "H_3(d)" << lcbopii.H_3(lcbopii.d) << std::endl;*/

   return 1;
}
