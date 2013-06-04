#include <iostream>
#include <Eigen/Dense>
#include <list>
#include <cmath>
#include "atom.h"
#include "lcbopii.h"
#include "printers.h"

using namespace  simul;


class EbSum: public AtomWalker
{
private:
	int n;
	double eb_total;
	LCBOPII lcbopii;

public:
	EbSum():
		n(0),
		eb_total(0.0)
	{
		//lcbopii = new LCBOPII();
	}

	~EbSum()
	{
	//	delete lcbopii;
	}

	double get_eb()
	{
		return eb_total/n;
	}

protected:
	virtual void bound_found(Atom *u, Atom *v)
	{
		std::cout << "Atom: " << u->get_id() << " --- " << v->get_id() << std::endl;
		//eb_total += lcbopii->Eb_ij(u,v);
	}
	virtual void atom_found(Atom *u)
	{
		std::cout << "Odwiedzilem Atom: " << u->get_id() << std::endl;
		Atom *v;
		for(Atom::bond_type::iterator it = u->get_bonds().begin();
				it != u->get_bonds().end(); it++)
		{
			v = *it;
			eb_total += lcbopii.Eb_ij(u,v);
		}
		n++;
	}
};

void create_di(std::list<Atom *> * atoms, double r)
{

	Atom * i = new Atom(0.0, 0, 0.0);
	Atom * j = new Atom(r, 0.0, 0.0);
	i->addBond(j);
	atoms->push_back(i);
	atoms->push_back(j);
}

void create_chi(std::list<Atom *> * atoms, double r)
{
	//double omega = 2.0*M_PI/3.0;

	Atom * i = new Atom(0.0, 0, 0.0);
	/*
	Atom * k1 = new Atom(r*std::cos(omega), r*std::sin(omega), 0.0);
	Atom * k2 = new Atom(r*std::cos(2*omega), r*std::sin(2*omega), 0.0);
	*/


	Atom * j = new Atom(r, 0.0, 0.0);

	Atom * k = new Atom(2*r, 0.0, 0.0);

	/*
	Atom * l1 = new Atom(r*std::cos(omega + M_PI) + r, r*std::sin(omega + M_PI), 0.0);
	Atom * l2 = new Atom(r*std::cos(2*omega + M_PI) + r, r*std::sin(2*omega + M_PI), 0.0);
	*/

	i->addBond(j);
	j->addBond(k);

	atoms->push_back(i);
	atoms->push_back(j);
	atoms->push_back(k);
}

void create_chi_z(std::list<Atom *> * atoms, double r)
{
	//double omega = 2.0*M_PI/3.0;

	Atom * i = new Atom(0.0, 0, 0.0);
	/*
	Atom * k1 = new Atom(r*std::cos(omega), r*std::sin(omega), 0.0);
	Atom * k2 = new Atom(r*std::cos(2*omega), r*std::sin(2*omega), 0.0);
	*/


	Atom * j = new Atom(r*std::cos(M_PI/3), r*std::sin(M_PI/3), 0.0);

	Atom * k = new Atom(r*std::cos(M_PI/3), -r*std::sin(M_PI/3), 0.0);

	/*
	Atom * l1 = new Atom(r*std::cos(omega + M_PI) + r, r*std::sin(omega + M_PI), 0.0);
	Atom * l2 = new Atom(r*std::cos(2*omega + M_PI) + r, r*std::sin(2*omega + M_PI), 0.0);
	*/

	i->addBond(j);
	i->addBond(k);

	atoms->push_back(i);
	atoms->push_back(j);
	atoms->push_back(k);
}

void create_tb(std::list<Atom *> * atoms, double r)
{

	Atom * i = new Atom(0.0, 0, 0.0);
	Atom * j = new Atom(r, 0.0, 0.0);

	Atom * k = new Atom(-r, 0.0, 0.0);
	Atom * l = new Atom(2*r, 0.0, 0.0);
	i->addBond(j);
	i->addBond(k);
	j->addBond(l);
	atoms->push_back(i);
	atoms->push_back(j);
	atoms->push_back(k);
	atoms->push_back(l);

	Atom * H1k = new Atom(-2*r, 0.0, 0.0);
	Atom * H2k = new Atom(-2*r, 0.0, 0.0);
	Atom * H3k = new Atom(-2*r, 0.0, 0.0);
	k->addBond(H1k);
	k->addBond(H2k);
	k->addBond(H3k);

	Atom * H1l = new Atom(3*r, 0.0, 0.0);
	Atom * H2l = new Atom(3*r, 0.0, 0.0);
	Atom * H3l = new Atom(3*r, 0.0, 0.0);
	l->addBond(H1l);
	l->addBond(H2l);
	l->addBond(H3l);

	atoms->push_back(H1k);
	atoms->push_back(H2k);
	atoms->push_back(H3k);
	atoms->push_back(H1l);
	atoms->push_back(H2l);
	atoms->push_back(H3l);
}

void create_db(std::list<Atom *> * atoms, double r)
{
	double omega = 2.0*M_PI/3.0;

	Atom * i = new Atom(0.0, 0, 0.0);
	Atom * j = new Atom(r, 0.0, 0.0);

	Atom * k1 = new Atom(r*std::cos(omega), r*std::sin(omega), 0.0);
	Atom * k2 = new Atom(r*std::cos(2*omega), r*std::sin(2*omega), 0.0);

	Atom * l1 = new Atom(r*std::cos(omega + M_PI) + r, r*std::sin(omega + M_PI), 0.0);
	Atom * l2 = new Atom(r*std::cos(2*omega + M_PI) + r, r*std::sin(2*omega + M_PI), 0.0);


	i->addBond(j);
	i->addBond(k1);
	i->addBond(k2);

	j->addBond(l1);
	j->addBond(l2);

	atoms->push_back(i);
	atoms->push_back(j);
	atoms->push_back(k1);
	atoms->push_back(k2);
	atoms->push_back(l1);
	atoms->push_back(l2);

	// H atoms
	Atom * H1k1 = new Atom(2*r*std::cos(omega), 2*r*std::sin(omega), 0.0);
	k1->addBond(H1k1);
	Atom * H2k1 = new Atom(2*r*std::cos(omega), 2*r*std::sin(omega), 0.0);
	k1->addBond(H2k1);
	Atom * H3k1 = new Atom(2*r*std::cos(omega), 2*r*std::sin(omega), 0.0);
	k1->addBond(H3k1);

	Atom * H1k2 = new Atom(2*r*std::cos(2*omega), 2*r*std::sin(2*omega), 0.0);
	k2->addBond(H1k2);
	Atom * H2k2 = new Atom(2*r*std::cos(2*omega), 2*r*std::sin(2*omega), 0.0);
	k2->addBond(H2k2);
	Atom * H3k2 = new Atom(2*r*std::cos(2*omega), 2*r*std::sin(2*omega), 0.0);
	k2->addBond(H3k2);

	Atom * H1l1 = new Atom(2*r*std::cos(omega + M_PI) + r, 2*r*std::sin(omega + M_PI), 0.0);
	l1->addBond(H1l1);
	Atom * H2l1 = new Atom(2*r*std::cos(omega + M_PI) + r, 2*r*std::sin(omega + M_PI), 0.0);
	l1->addBond(H2l1);
	Atom * H3l1 = new Atom(2*r*std::cos(omega + M_PI) + r, 2*r*std::sin(omega + M_PI), 0.0);
	l1->addBond(H3l1);

	Atom * H1l2 = new Atom(2*r*std::cos(2*omega + M_PI) + r, 2*r*std::sin(2*omega + M_PI), 0.0);
	l2->addBond(H1l2);
	Atom * H2l2 = new Atom(2*r*std::cos(2*omega + M_PI) + r, 2*r*std::sin(2*omega + M_PI), 0.0);
	l2->addBond(H2l2);
	Atom * H3l2 = new Atom(2*r*std::cos(2*omega + M_PI) + r, 2*r*std::sin(2*omega + M_PI), 0.0);
	l2->addBond(H3l2);

	atoms->push_back(H1k1);
	atoms->push_back(H2k1);
	atoms->push_back(H3k1);
	atoms->push_back(H1k2);
	atoms->push_back(H2k2);
	atoms->push_back(H3k2);
	atoms->push_back(H1l1);
	atoms->push_back(H2l1);
	atoms->push_back(H3l1);
	atoms->push_back(H1l2);
	atoms->push_back(H2l2);
	atoms->push_back(H3l2);
}

void create_db_rot(std::list<Atom *> * atoms, double r)
{
	double omega = 2.0*M_PI/3.0;

	Atom * i = new Atom(0.0, 0, 0.0);
	Atom * k1 = new Atom(r*std::cos(omega), r*std::sin(omega), 0.0);
	Atom * k2 = new Atom(r*std::cos(2*omega), r*std::sin(2*omega), 0.0);


	Atom * j = new Atom(r, 0.0, 0.0);
	Atom * l1 = new Atom(r*std::cos(omega + M_PI) + r, r*std::sin(omega + M_PI), 0.0);
	Atom * l2 = new Atom(r*std::cos(2*omega + M_PI) + r, r*std::sin(2*omega + M_PI), 0.0);

	Eigen::Vector3d ri = i->r - j->r;
	Eigen::Vector3d k1_org, k2_org;

	k1_org = k1->r;
	k2_org = k2->r;
	k1->r = k1_org;
	k2->r = k2_org;
	Eigen::AngleAxisd rot_ri(M_PI/2, ri.normalized());

	Eigen::Vector3d r_k1_rot = rot_ri * (k1->r - i->r);
	k1->r += (r_k1_rot - k1->r);

	Eigen::Vector3d r_k2_rot = rot_ri * (k2->r - i->r);
	k2->r += (r_k2_rot - k2->r);

	i->addBond(j);
	i->addBond(k1);
	i->addBond(k2);

	j->addBond(l1);
	j->addBond(l2);

	atoms->push_back(i);
	atoms->push_back(j);
	atoms->push_back(k1);
	atoms->push_back(k2);
	atoms->push_back(l1);
	atoms->push_back(l2);
}

void free_atoms(std::list<Atom *> * atoms)
{
	//for(std::list<Atom *>::iterator it = atoms->begin(); it != atoms->end(); ++it)
	while(atoms->size())
	{
		Atom * i = atoms->back();
		delete i;
		atoms->pop_back();
	}
}

double compute_Eb(std::list<Atom *> * atoms, int stop = -1)
{
	LCBOPII lcbopii;
	Atom * n_i;
	Atom * n_j;
	double eb, eb_total;
	int n = atoms->size();
	Atom::bond_type::iterator jt;
	Atom::bond_type tmp;
	eb_total = 0.0;
	if(stop <0)
		stop = n;

	for(std::list<Atom *>::iterator it = atoms->begin(); it != atoms->end(); ++it)
	{
		if(stop == 0) break;
		n_i = *it;
		tmp = n_i->get_bonds();
		jt = tmp.begin();
		std::cout << "i = " << n_i->get_id() << " " << n_i->get_bonds().size() << std::endl;
		while(jt != tmp.end())
		{
			n_j = *jt;
			std::cout << "j = " << n_i->get_id() << std::endl;
			if(n_j->get_id() == n_i->get_id()+1)
			{
			eb = lcbopii.Eb_ij(n_i,n_j);
			eb_total += eb;
			std::cout << n_j->r << std::endl;
			std::cout << "Eb_i = " << eb << std::endl;
			}
			//std::cout << std::endl;
			++jt;
		}
		stop--;
	}

	return eb_total;
}

int main()
{
	std::list<Atom *> atoms;

	for(double r = 0.9; r<1.65; r+=0.005)
	{
		create_tb(&atoms, r);
		std::cout << r << " LLL " << compute_Eb(&atoms,2) << std::endl << std::endl;
		free_atoms(&atoms);
	}
}
