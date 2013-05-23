#include <iostream>
#include <gmtl/gmtl.h>
#include "atom.h"

namespace simul
{
	int Atom::id_cnt = 0;

	Atom::Atom(double x, double y, double z):
		id(id_cnt++)
	{
		this->r = gmtl::Vec3d(x, y, z);
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
}
