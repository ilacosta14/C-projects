Uint CG(Sfield out, void (*op)(Sfield, Sfield), Sfield in, int size) {
    // Maximum number of iterations is defined in the code as CG_max_iter
    Scalar this_tol = CG_tolerance // Tolerance for convergence. Johannes used a somewhat more complicated definition. 
                                    // but I think that for the moment I can keep it like this

    // r is the residual: r = b - A*x
    Sfield r= new_Sfield(), p= new_Sfield(), Ap= new_Sfield();
    
    // for the moment, i choose to set x=out=0 initially
    set_Sfield2zero ( out );
    assign_Sfield ( r, in );

    // Initialize p with r for the first iteration
    
    assign_Sfield ( p, r );

    Scalar r_sq = dot_Sfield(r, r);
    Scalar alpha, beta, residual;
    Uint iter;

    //Start the CG iterations 
    for (iter = 1; iter <=  CG_max_iter; iter++) {
        // Compute A*p
        op(Ap, p);

        // Compute alpha = (r·r) / (p·Ap)
        alpha =  r_sq/ dot_Sfield(p, Ap);

        // Update the solution: x = x + alpha*p
        // this function takes  the field p, multiplies it by the factor alpha and adds it to the value of out
        inplace_Sfield_fma ( out, alpha, p );


        // Update the residual: r = r - alpha*Ap
        inplace_Sfield_fma ( r, -alpha, Ap );

        // Compute the new residual norm
        residual = sqrt(dot_Sfield(r, r));

        // Check for convergence
        if (residual < tolerance) {
            break;
        }

        // Compute beta = (new_r·new_r) / (r·r)
        Scalar old_r_sq=r_sq;
        r_sq = dot_Sfield(r, r)
        Scalar beta = r_sq / old_r_sq;

        // Update p: p = r + beta*p
        assign_Sfield_fma ( p, r, beta, p );

    }

    /* free temporary fields */
    destroy_Sfield ( &p );
    destroy_Sfield ( &Ap );
    destroy_Sfield ( &r );

    return iter; // Return the number of iterations performed
}
