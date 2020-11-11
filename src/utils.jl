function structure_err!(s1, s2)
    d1 = fetch_data(s1, :d)
    d2 = fetch_data(s2, :d)
    err = sum( map((x1,x2)->abs(x1-x2), d1, d2) ) / length(d1)
    # println("err = ",err)
    return err    
end