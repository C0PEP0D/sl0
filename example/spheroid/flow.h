#include <iostream>

// Simple shear flow
template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class Flow {
    public:
        Flow() : gamma(1.0) {}
    public:
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const;
        TypeVector getVorticity(const TypeRef<const TypeVector>& x, const double& t) const;
        TypeMatrix getStrain(const TypeRef<const TypeVector>& x, const double& t) const;
    public:
        double gamma;
};

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector Flow<TypeVector, TypeMatrix, TypeRef>::getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector u; u.fill(0.0);
    u(0) = gamma * x(1);
    return u;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeVector Flow<TypeVector, TypeMatrix, TypeRef>::getVorticity(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeVector o; o.fill(0.0);
    o(2) = -0.5 * gamma;
    return o;
}

template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
TypeMatrix Flow<TypeVector, TypeMatrix, TypeRef>::getStrain(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix strain; strain.fill(0.0);
    strain(0,1) = 0.5 * gamma;
    strain(1,0) = 0.5 * gamma;
    return strain;
}
