#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include "printers.h"
#include "atom.h"

const std::string SdfPrinter::INDENT = "  ";

std::string SdfPrinter::dump_atoms()
{
	std::stringstream atoms_ss;
	atoms_ss.precision(5);
	while (!atoms.empty())
	{
		atoms_ss
				 << std::setw(10)<< atoms.top()->r[0] << " "
				 << std::setw(10)<< atoms.top()->r[1] << " "
				 << std::setw(10)<< atoms.top()->r[2]
				 << " C ";
		for(int i=0; i<12; i++) atoms_ss << " 0 ";
		atoms_ss << std::endl;
	    atoms.pop();
	}
	return atoms_ss.str();
}
std::string SdfPrinter::dump_bonds()
{
	std::stringstream bonds_ss;
	std::pair<int,int> bond;
	while (!bonds.empty())
	{
		bond = bonds.top();
		bonds_ss
				 << std::setw(3) <<  bond.first << std::setw(3)  << bond.second
				 << "  1  0  0  0  0" << std::endl;
	    bonds.pop();
	}
	return bonds_ss.str();
}
std::string SdfPrinter::get_sdf(Atom * start)
{
	std::stringstream sdf_ss;
	walk_BFS(start);
	sdf_ss << "SDFPRINTER FOR LCBOPII" << std::endl;
	sdf_ss << " AUTHOR PAWEL TOMASIEWICZ" << std::endl;
	sdf_ss << std::endl;

	sdf_ss << std::setw(3) << atoms.size()  << std::setw(3) << " " << bonds.size()
			<< "  0  1  0  0  0  0  0  1 V2000" << std::endl;

	sdf_ss << dump_atoms();
	sdf_ss << dump_bonds();

	sdf_ss << "M END" << std::endl;

	return sdf_ss.str();
}
