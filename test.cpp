#include <gmtl/gmtl.h>
#include <gmtl/Matrix.h>
#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/shared_ptr.hpp>
#include <iostream>
#include <deque>
#include "atom.h"
#include "lcbopii.h"

using namespace  simul;
int main()
{
   //Atom::bond_type atoms, atoms2;
   std::deque<Atom *> atoms;
   for(int i=0; i<10; i++)
   {
	   /*
	   atoms.push_front(new Atom());
	   Atom::bond_type::iterator tmp = atoms.begin();
	   tmp->ppppp = i;
	   */
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
   }

   LCBOPII lcbopii;// = new LCBOPII();
   //int i = 0;
/*
   std::cout << lcbopii.S_down(0.0,-3.0) << std::endl;
   std::cout << lcbopii.S_down(0.5, -3.0) << std::endl;
   std::cout << lcbopii.S_down(1.0,-3.0) << std::endl;
*/

   gmtl::Vec3f v1(1.0, 0,0);
   gmtl::Vec3f v2(0.0, 1.0,0);
   gmtl::Vec3f diff = v1 - v2;

   float length = gmtl::length(diff);
   std::cout << length << std::endl;
   return 1;
}
