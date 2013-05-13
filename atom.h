#ifndef __ATOM_H
#define __ATOM_H
#include <set>
#include <iostream>
#include <gmtl/gmtl.h>


namespace simul
{
	class Atom
	{
	public:
		typedef Atom self_type;
		/*TODO zamienić to poniżej na map o indeksach atom.id i kluczach *atom */
		typedef std::set<Atom *> bond_type;
		typedef gmtl::Vec3d position_type;

	private:
		bond_type bonds;

	public:
		int id; //TODO przerobić na static z autoinkrementacja
		position_type r;

		Atom(double x, double y, double z);
		void addBond(Atom * j, bool reverse = true);
		int removeBond(Atom * j, bool reverse = true);

		bond_type get_bonds() const;
		int get_id() const;
	};
}
#endif
