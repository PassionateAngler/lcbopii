#ifndef __ATOM_H
#define __ATOM_H
#include <set>
#include <iostream>
#include <Eigen/Dense>


namespace simul
{
	class Atom
	{
	public:
		typedef Atom self_type;
		/*TODO zamienić to poniżej na map o indeksach atom.id i kluczach *atom */
		typedef std::set<Atom *> bond_type;
		typedef Eigen::Vector3d position_type;

	private:
		bond_type bonds;
		static int id_cnt;
		int id;

	public:
		position_type r;
		bool visited_switch;

		Atom(double x, double y, double z);
		~Atom()
		{
			//std::cout << id << " destroyed" << std::endl;
		}
		void addBond(Atom * j, bool reverse = true);
		int removeBond(Atom * j, bool reverse = true);

		bond_type get_bonds() const;
		int get_id() const;
	};

	class AtomWalker
	{
	public:
		/**
		 * Breadth-first search of atom structure graph
		 * @start - atom to start from
		 */
		void walk_BFS(Atom * start);
	protected:
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
