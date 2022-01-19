
template<typename TypeVector, template<typename...> typename TypeRef>
class Flow {
    public:
        TypeVector getVelocity(const TypeRef<const TypeVector>& x, const double& t) const {
            return x; // Exponential derivative
        }
};

