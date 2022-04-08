
class Flow {
    public:
        template<typename TypeVector>
        TypeVector getVelocity(const TypeVector& x, const double& t) const {
            return x; // Exponential derivative
        }
};

