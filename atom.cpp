#include <iostream>
#include <deque>
#include <Eigen/Dense>
#include "atom.h"

namespace simul
{
	int Atom::id_cnt = 0;

	Atom::Atom(double x, double y, double z):
		id(id_cnt++),
		visited_switch(false)
	{
		this->r = Eigen::Vector3d(x, y, z);
	}

	int Atom::removeBond(Atom *j, bool reverse)
	{
		Atom::bond_type::iterator it_j = this->bonds.find(j);
		if(it_j == this->bonds.end())
		{
			std::cerr << "Atom " << this->id
					  << " does not have neighbour " << j->id << std::endl;
			return -1;
		}

		this->bonds.erase(it_j);
		if(reverse)
		{
			j->removeBond(this, false);
		}
		return 0;
	}

	void Atom::addBond(Atom *j, bool reverse)
	{
		this->bonds.insert(j);
		if(reverse)
		{
			j->addBond(this, false);
		}
	}

	Atom::bond_type Atom::get_bonds() const
	{
		return bonds;
	}

	int Atom::get_id() const
	{
		return id;
	}

	void AtomWalker::walk_BFS(Atom * start)
	{
		Atom * u;
		Atom * v;
		std::deque<Atom *> Q;

		const bool NON_VISITED = start->visited_switch;
		const bool VISITED = !start->visited_switch;
		Q.push_back(start);
		while(!Q.empty())
		{
			u = Q.back();
			Q.pop_back();
			for(Atom::bond_type::iterator it = u->get_bonds().begin();
				it != u->get_bonds().end(); it++)
			{
				v = *it;
				if(v->visited_switch == NON_VISITED)
				{
					bound_found(u,v);
					Q.push_back(v);
				}
			}
			atom_found(u);
			u->visited_switch = VISITED;
		}
	}
}
