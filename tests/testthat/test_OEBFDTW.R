
set.seed(1)
data<-simulate_data(n_obs=100)

dom=c(0,1)
basis_x<-create.bspline.basis(c(0,1),nbasis =30)
x_fd<-smooth.basis(data$grid_i[[1]]/max(data$grid_i[[1]]),data$x_err[[1]],basis_x)$fd
template_fd<-smooth.basis(data$grid_template,data$template,basis_x)$fd
der_x_fd<-deriv.fd(x_fd,1)
der_template_fd=deriv.fd(template_fd,1)



mod<-OEBFDTW(x_fd,template_fd,der_x_fd ,der_template_fd,get_fd = "x_reg")
