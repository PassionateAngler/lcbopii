#include <iostream>
#include <gmtl/gmtl.h>
#include "atom.h"

namespace simul
{
	Atom::Atom(float x, float y, float z):
			id(-1)
	{
		this->r = gmtl::Vec3f(x, y, z);
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
