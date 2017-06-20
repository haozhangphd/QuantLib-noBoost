#ifndef quantlib_variategenerator_hpp
#define quantlib_variategenerator_hpp

template<class Engine, class Distribution>
class variate_generator {
public:
    using result_type = typename Distribution::result_type;

    variate_generator(Engine e, Distribution d)
            : eng{e}, dist{d} {}

    result_type operator()() { return dist(eng); }

    template<class T>
    result_type operator()(const T &value) { return dist(eng, value); }

    const Engine &engine() const { return eng; }

    const Distribution &distribution() const { return dist; }


private:
    Engine eng;
    Distribution dist;
};


#endif
