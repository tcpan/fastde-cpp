#include "utils_data.tpp"


// ------- explicit instantiation

template cpp11::writable::doubles export_vec_to_rvec(double * pv, size_t const & count);
template cpp11::writable::doubles export_vec_to_rvec(const double * pv, size_t const & count);
template cpp11::writable::logicals export_vec_to_rvec(bool * pv, size_t const & count);
template cpp11::writable::logicals export_vec_to_rvec(const bool * pv, size_t const & count);

// template cpp11::writable::doubles export_vec_to_rvec(std::vector<double> const & pv);
// template cpp11::writable::logicals export_vec_to_rvec(std::vector<bool> const & pv);



template cpp11::writable::logicals_matrix<cpp11::by_column> export_vec_to_r_matrix(
    std::vector<bool> const & pv, size_t const & nrow, size_t const & ncol
);
template cpp11::writable::logicals_matrix<cpp11::by_column> export_vec_to_r_matrix(
    bool * pv, size_t const & nrow, size_t const & ncol
);
template cpp11::writable::logicals_matrix<cpp11::by_column> export_vec_to_r_matrix(
    const bool * pv, size_t const & nrow, size_t const & ncol
);


template cpp11::writable::doubles_matrix<cpp11::by_column> export_vec_to_r_matrix(
    std::vector<double> const & pv, size_t const & nrow, size_t const & ncol
);
template cpp11::writable::doubles_matrix<cpp11::by_column> export_vec_to_r_matrix(
    double * pv, size_t const & nrow, size_t const & ncol
);
template cpp11::writable::doubles_matrix<cpp11::by_column> export_vec_to_r_matrix(
    const double * pv, size_t const & nrow, size_t const & ncol
);


template cpp11::writable::list export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels
);

template cpp11::writable::list export_fc_to_r_matrix(
    double * fc, std::string const & fcname,
    double * p1, std::string const & p1name,
    double * p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels
);
template cpp11::writable::list export_fc_to_r_matrix(
    const double * fc, std::string const & fcname,
    const double * p1, std::string const & p1name,
    const double * p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels
);


template cpp11::writable::data_frame export_vec_to_r_dataframe(
    std::vector<double> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
template cpp11::writable::data_frame export_vec_to_r_dataframe(
    double * pv, std::string const & name, size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
template cpp11::writable::data_frame export_vec_to_r_dataframe(
    const double * pv, std::string const & name, size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);


template cpp11::writable::data_frame export_fc_to_r_dataframe(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
template cpp11::writable::data_frame export_fc_to_r_dataframe(
    double * fc, std::string const & fcname,
    double * p1, std::string const & p1name,
    double * p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);
template cpp11::writable::data_frame export_fc_to_r_dataframe(
    const double * fc, std::string const & fcname,
    const double * p1, std::string const & p1name,
    const double * p2, std::string const & p2name,
    size_t const & count,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    cpp11::strings const & features
);

