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
		size_t n = r.size();
		array1d<T> x(n);
		
		for (size_t i=1; i<n; i++) {
				T w = a(i) / b(i-1);
				b(i) -= w * c(i-1);
				r(i) -= w * r(i-1);
		}
		
		x(n-1) = r(n-1) / b(n-1);
		
		for (size_t i=n-2; i>0; i--) {
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
}


class finite_difference_AF {
		public:

			using scalar = finite_differences::complex;
			using field = finite_differences::array_1D;

		private:

			field up,rfp,rap;

		public:

			field u,rf,ra;

			void step(){
				size_t n = u.size()-2;
				
				finite_differences::array_1D A(n), B(n), R(n);

				A = - ra.segment(n, 1) / 2;
				B = 1 + ra.segment(n,1) - rf.segment(n,1);
				R = (up.segment(n,2) + up.segment(n,0)) * ra.segment(n,1) / 2
				  + up.segment(n,1) * (1 + rfp.segment(n,1) - rap.segment(n,1));

				R(0)   += u(0) * ra(0)/2.;
				R(n-1) += u(n+1) * ra(n+1)/2.;

				u.segment(n, 1) = algebra::tridiagonal(A,B,A,R);
			};
			
			void update(){
				std::swap(ra, rap);
				std::swap(rf, rfp);
				std::swap(u, up);
			};

			void resize(size_t N){
				rf.resize(N);
				ra.resize(N);
				rfp.resize(N);
				rap.resize(N);
				u.resize(N);
				up.resize(N);
			}
	};

	class finite_difference_ACF{
		public:

			using scalar = finite_differences::complex;
			using field = finite_differences::array_2D;

		private:

			field up,rfp,rap,rcp;


			template<typename field>
			void ACF_step(const field &ra, const field &rc, const field &rf, const field &rap, const field &rcp, const field &rfp, field &u, const field &up){

				size_t nx = u.cols()-2;
				size_t ny = u.rows()-2;

				// TODO: parallel
				
				finite_differences::array_1D A(ny), B(ny), R(ny);
				for(size_t i=1; i<nx+1; i++){
					
					const auto ra_i = ra.col(i);
//					const auto rc_i = rc.col(i);
					const auto rf_i = rf.col(i);
//					const auto rap_i = rap.col(i);
					const auto rcp_i = rcp.col(i);
					const auto rfp_i = rfp.col(i);
					auto u_i = u.col(i);
					const auto up_i = up.col(i);
					
					A = -ra_i.segment(ny, 1);
					B = 1 + 2 * ra_i.segment(ny, 1) - rf_i.segment(ny, 1);
					R = (up_i.segment(ny, 2) + up_i.segment(ny, 0)) * rcp_i.segment(ny,1)
					  + up_i.segment(ny,1) * (1 + rfp_i.segment(ny,1) - 2 * rcp_i.segment(ny,1));

					R(0) += u_i(0) * ra_i(0);
					R(ny-1) += u_i(ny+1) * ra_i(ny+1);

					u_i.segment(ny,1) = algebra::tridiagonal(A, B, A, R);
				}
			}

		public:

			field u, rf, ra, rc;
			unsigned thread_count = 1;

			void step_1(){
				ACF_step(ra, rc, rf, rap, rcp, rfp, u, up);
			};
			void step_2(){
				auto u_transposed = u.transpose();
				ACF_step(rc.transpose(), ra.transpose(), rf.transpose(), rcp.transpose(), rap.transpose(), rfp.transpose(), u_transposed, up.transpose());
			};
			void update(){
				std::swap(rf, rfp);
				std::swap(ra, rap);
				std::swap(rc, rcp);
				std::swap(u, up);
			};

			void resize(size_t Nx, size_t Ny){
				ra.resize(Nx,Ny);
				rap.resize(Nx,Ny);
				rc.resize(Nx,Ny);
				rcp.resize(Nx,Ny);
				rf.resize(Nx,Ny);
				rfp.resize(Nx,Ny);
				u.resize(Nx,Ny);
				up.resize(Nx,Ny);
			};
	};


	class finite_difference_A0F{
		public:

			using scalar = finite_differences::complex;
			using field = finite_differences::array_2D;

		private:

			field up, rfp, rap;

		public:

			field u, rf, ra;
			
			unsigned thread_count = 1;

			void step(){
				unsigned nx = u.rows()-2;
				unsigned ny = u.cols();
				
				// TODO: parallel
				finite_differences::array_1D A(nx), B(nx), R(nx);
				for(size_t i=0; i<ny; i++){
					
						const auto ra_i = ra.col(i);
						const auto rf_i = rf.col(i);
						const auto up_i = up.col(i);
						const auto rap_i = rap.col(i);
						const auto rfp_i = rfp.col(i);
						auto u_i = u.col(i);

						A = -ra_i.segment(nx,1) / 2.;
						B = 1 + ra_i.segment(nx, 1) - rf_i.segment(nx, 1);
						R = (up_i.segment(nx, 2) + up_i.segment(nx,0)) * rap_i.segment(nx,1)
						    + up_i.segment(nx, 1) * (1 + rfp_i.segment(nx, 1) - rap_i.segment(nx, 1));

						R(0) += u_i(0) * ra_i(0) / 2.;
						R(nx-1) += u_i(nx+1) * ra_i(nx+1) / 2.;

						u_i.segment(nx,1) = algebra::tridiagonal(A, B, A, R);
				}
			};
			
			
			void update(){
				std::swap(ra, rap);
				std::swap(rf, rfp);
				std::swap(u, up);
			};


			void resize(size_t Nx, size_t Ny){
				rf.resize(Nx,Ny);
				ra.resize(Nx,Ny);
				rfp.resize(Nx,Ny);
				rap.resize(Nx,Ny);
				u.resize(Nx,Ny);
				up.resize(Nx,Ny);
			};
	};


	class finite_difference_ABC{
		public:

			using scalar = finite_differences::complex;
			using field = finite_differences::array_2D;

		private:

			field up, rap, rbp, rcp, rzp; 

		public:
			// alpha dz u = A dr^2 u + B d_r u + C u 
			// ra = A/2Dr^2, rb = B/4Dr, rc = C/2, rz = alpha/Dz
			field u, ra, rb, rc, rz; 
			unsigned thread_count = 1;

			void step(){
				size_t nx = u.rows()-2;
				size_t ny = u.cols();
				
				finite_differences::array_1D A(nx), B(nx), C(nx), R(nx);

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

						A = rb_i.segment(nx, 1) - ra_i.segment(nx,1);
						B = 2 * ra_i.segment(nx, 1) + rz_i.segment(nx, 1) - rc_i.segment(nx, 1);
						C = - rb_i.segment(nx, 1) - ra_i.segment(nx, 1);

						R = (rap_i.segment(nx,1) + rbp_i.segment(nx,1)) * up_i.segment(nx,2)
						  + (rap_i.segment(nx,1) - rbp_i.segment(nx,1)) * up_i.segment(nx,0)
						  + (rcp_i.segment(nx,1) + rzp_i.segment(nx,1) - 2 * rap_i.segment(nx,1)) * up_i.segment(nx,1); 

						R(0)    += (-rb_i(1) + ra_i(1)) * u_i(0);
						R(nx-1) += ( rb_i(nx) + ra_i(nx)) * u_i(nx+1);

						u_i.segment(nx, 1) = algebra::tridiagonal(A, B, C, R);
				}

			};
			
			
			void update(){
				std::swap(ra, rap);
				std::swap(rb, rbp);
				std::swap(rc, rcp);
				std::swap(rz, rzp);
				std::swap(u, up);
			};
			

		void resize(size_t Nx,size_t Ny=1){
				ra.resize(Nx,Ny);
				rb.resize(Nx,Ny);
				rc.resize(Nx,Ny);
				rz.resize(Nx,Ny);
				rap.resize(Nx,Ny);
				rbp.resize(Nx,Ny);
				rcp.resize(Nx,Ny);
				rzp.resize(Nx,Ny);
				u.resize(Nx,Ny);
				up.resize(Nx,Ny);
		};
};
