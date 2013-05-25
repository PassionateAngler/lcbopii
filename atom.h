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
		static int id_cnt;
		int id;
		typedef enum{WHITE, GREY, BLACK} walk_type;
		unsigned long long walk_status;

	public:
		position_type r;
		bool visited_switch;

		Atom(double x, double y, double z);
		void addBond(Atom * j, bool reverse = true);
		int removeBond(Atom * j, bool reverse = true);

		bond_type get_bonds() const;
		int get_id() const;
	};

	class AtomWalker
	{
	protected:
		/**
		 * Breadth-first search of atom structure graph
		 * @start - atom to start from
		 */
		void walk_BFS(Atom * start);
		virtual void bound_found(Atom *u, Atom *v)
		{
			std::cout << u->get_id() << " --- " << v->get_id() << std::endl;
		};
		virtual void atom_found(Atom *u)
		{
			std::cout << u->get_id() << " : " << u->r << std::endl;
		};
		virtual ~AtomWalker(){}
	};
}
#endif
