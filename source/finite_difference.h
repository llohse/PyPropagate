#pragma once

#include <complex>
#include <Eigen/Dense>

namespace algebra{

	template <typename T>
		using array1d = Eigen::Array<T, Eigen::Dynamic, 1>;

	/* Solve the tridiagonal matrix equation
	 * a_i x_{i-1} + b_i x_i  + c_i x_{i+1} = r_i
	 * 
	 * Returns x
	 * 
	 * Modifies the input coefficients a, b, c, r!
	 */
	template <typename T>
		array1d<T> tridiagonal(array1d<T> &a, array1d<T> &b, array1d<T> &c, array1d<T> &r) {
			int n = r.size();
			array1d<T> x(n,1);

			for (int i=1; i<n; i++) {
				T w = a(i) / b(i-1);
				b(i) -= w * c(i-1);
				r(i) -= w * r(i-1);
			}

			x(n-1) = r(n-1) / b(n-1);

			for (int i=n-2; i>=0; i--) {
				x(i) = (r(i) - c(i) * x(i+1)) / b(i);
			}

			return x;
		}
}

namespace finite_differences{

	using real = double;
	using complex = std::complex<real>;

	using scalar = complex;

	using array_1D   = Eigen::Array<scalar, Eigen::Dynamic, 1>;
	using array_2D   = Eigen::Array<scalar, Eigen::Dynamic, Eigen::Dynamic>;

	using shape_t = std::pair<size_t, size_t>;

	template <typename T>
		shape_t get_shape(const Eigen::EigenBase<T> &a) {
			return {a.rows(), a.cols()};
		}

	template <typename A>
		void set_if_compatible(Eigen::DenseBase<A> & dst, const Eigen::Ref<const A> & src) {
			if (get_shape(src) != get_shape(dst))
				throw std::invalid_argument("shape does not match");
			dst = src;
		}


	class finite_difference_AF {
		public:

			using scalar = complex;
			using field = array_1D;

		private:

			field up,rfp,rap;

		public:

			field u,rf,ra;


			finite_difference_AF  (size_t n) 
				: up(n,1), rfp(n,1), rap(n,1), u(n,1), rf(n,1), ra(n,1) {};

			void step(){
				size_t n = u.size()-2;

				field A = - ra.segment(1, n) / 2;
				field B = 1 + ra.segment(1, n) - rf.segment(1, n);
				field C = A;
				field R = (up.segment(2,n) + up.segment(0,n)) * ra.segment(1,n) / 2
					+ up.segment(1,n) * (1 + rfp.segment(1,n) - rap.segment(1,n));

				R(0)   += u(0) * ra(0)/2.;
				R(n-1) += u(n+1) * ra(n+1)/2.;

				u.segment(1,n) = algebra::tridiagonal(A,B,C,R);
			};

			void update(){
				std::swap(ra, rap);
				std::swap(rf, rfp);
				std::swap(u, up);
			};

			field & get_field() {
				return u;
			}


			void set_field(const Eigen::Ref<const field> v) {
				set_if_compatible(u, v);
			}

			void set_ra(const Eigen::Ref<const field> v) {
				set_if_compatible(ra, v);
			}

			void set_rf(const Eigen::Ref<const field> v) {
				set_if_compatible(rf, v);
			}

	};

	class finite_difference_ACF{
		public:

			using scalar = complex;
			using field = array_2D;

		private:

			field up,rfp,rap,rcp;


			template<typename field>
				static field ACF_step(const field & ra_a, 
						const field & rc_a, 
						const field & rf_a, 
						const field & rap_a, 
						const field & rcp_a, 
						const field & rfp_a, 
						const field & u_a, 
						const field & up_a){

					size_t nx = u_a.cols()-2;
					size_t ny = u_a.rows()-2;

					field res = u_a;

					// TODO: parallel

					array_1D A(ny), B(ny), C(ny), R(ny);
					for(size_t i=1; i<nx+1; i++){

						const auto ra_i = ra_a.col(i);
						const auto rf_i = rf_a.col(i);
						const auto rcp_i = rcp_a.col(i);
						const auto rfp_i = rfp_a.col(i);
						const auto u_i = u_a.col(i);

						A = -ra_i.segment(1,ny);
						B = 1 + 2 * ra_i.segment(1,ny) - rf_i.segment(1,ny);
						C = A;
						R = (up_a.col(i+1).segment(1,ny) + up_a.col(i-1).segment(1,ny)) * rcp_i.segment(1,ny)
							+  up_a.col(i).segment(1,ny) * (1 + rfp_i.segment(1,ny) - 2 * rcp_i.segment(1,ny));

						R(0) += u_i(0) * ra_i(0);
						R(ny-1) += u_i(ny+1) * ra_i(ny+1);

						res.col(i).segment(1,ny) = algebra::tridiagonal(A, B, C, R);
					}

					return res;
				}

		public:

			field u, rf, ra, rc;
			unsigned thread_count = 1;

			finite_difference_ACF(size_t ny, size_t nx) 
				: up(ny, nx), rfp(ny,nx), rap(ny,nx), rcp(ny,nx), u(ny,nx), rf(ny,nx), ra(ny,nx), rc(ny,nx) {};

			void step_1(){
				u = ACF_step(ra, rc, rf, rap, rcp, rfp, u, up);
			};

			void step_2(){
				//				auto u_transposed = u.transpose();
				u = ACF_step(rc.transpose(), ra.transpose(), rf.transpose(), rcp.transpose(), rap.transpose(), rfp.transpose(), u.transpose(), up.transpose()).transpose();
			};

			void update(){
				std::swap(rf, rfp);
				std::swap(ra, rap);
				std::swap(rc, rcp);
				std::swap(u, up);
			};

			field & get_field() {
				return u;
			}

			void set_field(const Eigen::Ref<const field> v) {
				set_if_compatible(u, v);
			}

			void set_ra(const Eigen::Ref<const field> v) {
				set_if_compatible(ra, v);
			}

			void set_rc(const Eigen::Ref<const field> v) {
				set_if_compatible(rc, v);
			}

			void set_rf(const Eigen::Ref<const field> v) {
				set_if_compatible(rf, v);
			}

	};


	class finite_difference_A0F{
		public:

			using scalar = complex;
			using field = array_2D;

		private:

			field up, rfp, rap;

		public:

			field u, rf, ra;

			unsigned thread_count = 1;

			finite_difference_A0F(size_t ny, size_t nx) 
				: up(ny,nx), rfp(ny,nx), rap(ny,nx), u(ny,nx), rf(ny,nx), ra(ny,nx) {}


			void step(){
				unsigned nx = u.cols();
				unsigned ny = u.rows()-2;

				// TODO: parallel
				array_1D A(ny), B(ny), C(ny), R(ny);
				for(size_t i=0; i<nx; i++){

					const auto ra_i = ra.col(i);
					const auto rf_i = rf.col(i);
					const auto up_i = up.col(i);
					const auto rap_i = rap.col(i);
					const auto rfp_i = rfp.col(i);
					auto u_i = u.col(i);

					A = -ra_i.segment(1,ny) / 2.;
					B = 1 + ra_i.segment(1,ny) - rf_i.segment(1,ny);
					C = A;
					R = (up_i.segment(2,ny) + up_i.segment(0,ny)) * rap_i.segment(1,ny)
						+ up_i.segment(1,ny) * (1 + rfp_i.segment(1,ny) - rap_i.segment(1,ny));

					R(0) += u_i(0) * ra_i(0) / 2.;
					R(ny-1) += u_i(ny+1) * ra_i(ny+1) / 2.;

					u_i.segment(1,ny) = algebra::tridiagonal(A, B, A, R);
				}
			}


			void update(){
				std::swap(ra, rap);
				std::swap(rf, rfp);
				std::swap(u, up);
			}

			field & get_field() {
				return u;
			}

			void set_field(const Eigen::Ref<const field> v) {
				set_if_compatible(u, v);
			}

			void set_ra(const Eigen::Ref<const field> v) {
				set_if_compatible(ra, v);
			}

			void set_rf(const Eigen::Ref<const field> v) {
				set_if_compatible(rf, v);
			}

	};


	class finite_difference_ABC{
		public:

			using scalar = complex;
			using field = array_2D;

		private:

			field up, rap, rbp, rcp, rzp; 

		public:
			field u, ra, rb, rc, rz; 
			// alpha dz u = A dr^2 u + B d_r u + C u 
			// ra = A/2Dr^2, rb = B/4Dr, rc = C/2, rz = alpha/Dz
			unsigned thread_count = 1;

			finite_difference_ABC(size_t ny, size_t nx)
				: up(ny,nx), rap(ny,nx), rbp(ny,nx), rcp(ny,nx), rzp(ny,nx), 
				u(ny,nx), ra(ny,nx), rb(ny,nx), rc(ny,nx), rz(ny,nx) {}

			void step(){
				size_t nx = u.rows()-2;
				size_t ny = u.cols();

				array_1D A(nx), B(nx), C(nx), R(nx);

				// TODO: parallel
				for (size_t i=0; i<ny; i++){

					const auto ra_i = ra.col(i);
					const auto rb_i = rb.col(i);
					const auto rc_i = rc.col(i);
					const auto rz_i = rz.col(i);

					const auto rap_i = rap.col(i);
					const auto rbp_i = rbp.col(i);
					const auto rcp_i = rcp.col(i);
					const auto rzp_i = rzp.col(i);

					const auto up_i = up.col(i);
					auto u_i = u.col(i);

					A = rb_i.segment(1,nx) - ra_i.segment(1,nx);
					B = 2 * ra_i.segment(1,nx) + rz_i.segment(1,nx) - rc_i.segment(1,nx);
					C = - rb_i.segment(1,nx) - ra_i.segment(1,nx);

					R = (rap_i.segment(1,nx) + rbp_i.segment(1,nx)) * up_i.segment(2,nx)
						+ (rap_i.segment(1,nx) - rbp_i.segment(1,nx)) * up_i.segment(0,nx)
						+ (rcp_i.segment(1,nx) + rzp_i.segment(1,nx) - 2 * rap_i.segment(1,nx)) * up_i.segment(1,nx); 

					R(0)    += (-rb_i(1) + ra_i(1)) * u_i(0);
					R(nx-1) += ( rb_i(nx) + ra_i(nx)) * u_i(nx+1);

					u_i.segment(1, nx) = algebra::tridiagonal(A, B, C, R);
				}

			}


			void update(){
				std::swap(ra, rap);
				std::swap(rb, rbp);
				std::swap(rc, rcp);
				std::swap(rz, rzp);
				std::swap(u, up);
			}

			field & get_field() {
				return u;
			}			

			void set_field(const Eigen::Ref<const field> v) {
				set_if_compatible(u, v);
			}

			void set_ra(const Eigen::Ref<const field> v) {
				set_if_compatible(ra, v);
			}

			void set_rb(const Eigen::Ref<const field> v) {
				set_if_compatible(rb, v);
			}

			void set_rc(const Eigen::Ref<const field> v) {
				set_if_compatible(rc, v);
			}

			void set_rz(const Eigen::Ref<const field> v) {
				set_if_compatible(rz, v);
			}
	};

}
