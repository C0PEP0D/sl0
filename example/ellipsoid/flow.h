#include <iostream>

// Simple shear flow
template<typename TypeVector, typename TypeMatrix, template<typename...> class TypeRef>
class Flow {
    public:
        Flow() : gamma(1.0) {}
    public:
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const;
        TypeMatrix getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const;
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
TypeMatrix Flow<TypeVector, TypeMatrix, TypeRef>::getVelocityGradients(const TypeRef<const TypeVector>& x, const double& t) const {
    TypeMatrix grads = TypeMatrix::Zero(); 
    grads(0,1) = gamma;
    return grads;
}
