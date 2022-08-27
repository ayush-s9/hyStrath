// TODO: Implement hllc flux function
//------------------------------------------------------------------------------
//  Harten-Lax-vanLeer with Contact flux function
//------------------------------------------------------------------------------
void HLLCFlux(
    const bool &downgradeSingleT,
    const bool &downgradeSingleTv,
    const bool &meshMoving,
    const scalar &mesh_phi,
    const vector &Sf,
    const scalar &cSf_pos,
    const scalar &cSf_neg,
    const scalar &rho_pos,
    const scalar &rho_neg,
    const PtrList<scalar> &pD_pos,
    const PtrList<scalar> &pD_neg,
    const vector &U_pos,
    const vector &U_neg,
    const vector &rhoU_pos,
    const vector &rhoU_neg,
    const scalar &p_pos,
    const scalar &p_neg,
    const PtrList<scalar> &evk_pos,
    const PtrList<scalar> &evk_neg,
    const scalar &e_pos,
    const scalar &e_neg,
    scalar &amaxSf,
    vector &aByU,
    scalar &rhoFlux,
    vector &rhoUFlux,
    PtrList<scalar> &rhoEvkFlux,
    scalar &rhoEFlux)
{
    scalar phiv_pos = U_pos & Sf;
    scalar phiv_neg = U_neg & Sf;

    //- Make fluxes relative to mesh-motion
    if (meshMoving)
    {
        phiv_pos -= mesh_phi;
        phiv_neg -= mesh_phi;
    }

    scalar u_pos = phiv_pos / mag(Sf);
    scalar u_neg = phiv_neg / mag(Sf);

    scalar c_pos = cSf_pos / mag(Sf); //- wave speed estimate +
    scalar c_neg = cSf_neg / mag(Sf); //- wave speed estimate -
    scalar cf = 0.5 * (c_pos + c_neg);

    scalar rhof = 0.5 * (rho_pos + rho_neg);

    // scalar pC = 0.5 * (p_pos + p_neg) + 0.5 * (u_neg - u_pos) * rhof * cf;
    scalar uC = 0.5 * (u_neg + u_pos) + 0.5 * (p_neg - p_pos) / (rhof * cf);

    // Can be changed for testing or can be made modular by implementing an additional option in the fvScheme dictionary
    scalar s_neg = min(u_neg - c_neg, u_pos - c_pos);
    scalar s_pos = max(u_neg + c_neg, u_pos + c_pos);

    if (0 <= s_neg)
    {
        rhoFlux = rho_neg * phiv_neg;
        rhoUFlux = rhoU_neg * phiv_neg + p_neg * Sf;
        rhoEFlux = phiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg);

        //- Make flux for pressure-work absolute
        if (meshMoving)
        {
            rhoEFlux += mesh_phi * p_neg;
        }

        if (downgradeSingleT)
        {
            rhoEvkFlux.clear();
        }
        else if (downgradeSingleTv)
        {
            rhoEvkFlux[0] = phiv_neg * rho_neg * evk_neg[0];
        }
        else
        {
            forAll(rhoEvkFlux, speciei)
            {
                rhoEvkFlux[speciei] = phiv_neg * (pD_neg[speciei] * evk_neg[speciei]);
            }
        }
    }

    if (u_pos <= 0)
    {
        rhoFlux = rho_pos * phiv_pos;
        rhoUFlux = rhoU_pos * phiv_pos + p_pos * Sf;
        rhoEFlux = phiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos);

        //- Make flux for pressure-work absolute
        if (meshMoving)
        {
            rhoEFlux += mesh_phi * p_pos;
        }

        if (downgradeSingleT)
        {
            rhoEvkFlux.clear();
        }
        else if (downgradeSingleTv)
        {
            rhoEvkFlux[0] = phiv_pos * rho_pos * evk_pos[0];
        }
        else
        {
            forAll(rhoEvkFlux, speciei)
            {
                rhoEvkFlux[speciei] = phiv_pos * (pD_pos[speciei] * evk_pos[speciei]);
            }
        }
    }

    if (s_neg <= 0 and 0 <= uC)
    {
        rhoFlux = rho_neg * phiv_neg + mag(Sf) * s_neg * (rho_neg * (s_neg - u_neg) / (s_neg - uC) - rho_neg);
        rhoUFlux = rhoU_neg * phiv_neg + p_neg * Sf + Sf * s_neg * (rho_neg * (s_neg - u_neg) / (s_neg - uC) - rho_neg) * uC;
        rhoEFlux = phiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + rho_neg * mag(Sf) * s_neg * (rho_neg * (s_neg - u_neg) / (s_neg - uC) - rho_neg) * (e_neg / rho_neg - (uC - u_neg) * (uC + p_neg / (rho_neg * (s_neg - u_neg))));

        //- Make flux for pressure-work absolute
        if (meshMoving)
        {
            rhoEFlux += mesh_phi * p_neg;
        }

        if (downgradeSingleT)
        {
            rhoEvkFlux.clear();
        }
        else if (downgradeSingleTv)
        {
            rhoEvkFlux[0] = phiv_neg * rho_neg * evk_neg[0] + rho_neg * mag(Sf) * s_neg * (evk_neg[0] * (s_neg - u_neg) / (s_neg - uC) - evk_neg[0]);
        }
        else
        {
            forAll(rhoEvkFlux, speciei)
            {
                rhoEvkFlux[speciei] = phiv_neg * (pD_neg[speciei] * evk_neg[speciei]) + pD_neg[speciei] * mag(Sf) * s_neg * (evk_neg[speciei] * (s_neg - u_neg) / (s_neg - uC) - evk_neg[speciei]);;
            }
        }
    }

    if (uC <= 0 and 0 <= s_pos)
    {
        rhoFlux = rho_pos * phiv_pos + mag(Sf) * s_pos * (rho_pos * (s_pos - u_pos) / (s_pos - uC) - rho_pos);
        rhoUFlux = rhoU_pos * phiv_pos + p_pos * Sf + Sf * s_pos * (rho_pos * (s_pos - u_pos) / (s_pos - uC) - rho_pos) * uC;
        rhoEFlux = phiv_pos * (rho_pos * (e_pos + 0.5 * magSqr(U_pos)) + p_pos) + rho_pos * mag(Sf) * s_pos * (rho_pos * (s_pos - u_pos) / (s_pos - uC) - rho_pos) * (e_pos / rho_pos - (uC - u_pos) * (uC + p_pos / (rho_pos * (s_pos - u_pos))));

        //- Make flux for pressure-work absolute
        if (meshMoving)
        {
            rhoEFlux += mesh_phi * p_pos;
        }

        if (downgradeSingleT)
        {
            rhoEvkFlux.clear();
        }
        else if (downgradeSingleTv)
        {
            rhoEvkFlux[0] = phiv_pos * rho_pos * evk_pos[0] + rho_pos * mag(Sf) * s_pos * (evk_pos[0] * (s_pos - u_pos) / (s_pos - uC) - evk_pos[0]);
        }
        else
        {
            forAll(rhoEvkFlux, speciei)
            {
                rhoEvkFlux[speciei] = phiv_pos * (pD_pos[speciei] * evk_pos[speciei]) + pD_pos[speciei] * mag(Sf) * s_pos * (evk_pos[speciei] * (s_pos - u_pos) / (s_pos - uC) - evk_pos[speciei]);
            }
        }
    }
}