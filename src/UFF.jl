module UFF

export bond_rest_length
export bond_energy
export angle_energy
export dihedral_energy
export vdw_energy
export electrostatic_energy

function bond_rest_length(atom1, atom2, bond_order)
    lambda = 0.13332
    n = bond_order
    r_BO = -lambda*(atom1.r1 + atom2.r1) * log(n)
    r_EN = atom1.r1*atom2.r1*(sqrt(atom1.Xi)-sqrt(atom2.Xi))^2/(atom1.Xi*atom1.r1+atom2.Xi*atom2.r1)
    r_IJ = atom1.r1 + atom2.r1 + r_BO + r_EN
    return r_IJ
end

function bond_energy(r, bond_order, atom1, atom2)
    r_IJ = bond_rest_length(atom1, atom2, bond_order)
    k_IJ = 664.12*atom1.Z1*atom2.Z1/(r_IJ^3)
    bond_energy = 1/2*k_IJ*(r - r_IJ)^2
    return bond_energy
end

function angle_energy(theta, atom1, atom2, atom3, bond_order_1, bond_order_2, bondtype)
    # Takes in theta in radians
    theta0 = atom2.theta0*pi/180
    r_IJ = bond_rest_length(atom1, atom2, bond_order_1)
    r_JK = bond_rest_length(atom2, atom3, bond_order_2)
    r_IK = sqrt(r_IJ^2 + r_JK^2 - 2*r_IJ*r_JK*cos(theta))
    beta = 664.12/r_IJ/r_JK
    K_IJK = beta*atom1.Z1*atom3.Z1/(r_IK^5)*r_IJ*r_JK*(3*r_IJ*r_JK*(1-cos(theta0)^2)-r_IK^2*cos(theta0))
    if bondtype == 0
        C_2 = 1/(4*sin(theta0)^2)
        C_1 = -4*C_2*cos(theta0)
        C_0 = C_2*(2*cos(theta0)^2+1)
        return K_IJK*(C_0+C_1*cos(theta)+C_2*cos(2*theta))
    end
    return K_IJK/bondtype^2*(1-cos(bondtype*theta))
end

function dihedral_energy(phi, atom2, atom3, atom2_is_sp3, atom3_is_sp3, n, phi0, central_bond_order)
    # take in phi and phi0 in radians
    if atom2_is_sp3 && atom3_is_sp3
        V = sqrt(atom2.Vi*atom3.Vi)
    else
        V = 5*sqrt(atom2.Uj*atom3.Uj)*(1+4.18*log(central_bond_order))
    end
    return 1/2*V*(1-cos(n*phi0)*cos(n*phi))
end

function vdw_energy(x, atom1, atom2)
    D_IJ = sqrt(atom1.D1*atom2.D1)
    x_IJ = sqrt(atom1.x1*atom2.x1)
    return D_IJ*(-2*(x_IJ/x)^6 + (x_IJ/x)^12)
end

function electrostatic_energy(r, q1, q2)
    return 332.0637*q1*q2/1/r
end

end
