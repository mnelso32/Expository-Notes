function B = tobinarytest(n)

	i = 1; 	
	q = floor(n/2); 
	r = n - 2*q;
	B(1) = r; 
	while 2 <= q     
		n = q;     
		i = i + 1;     
		q = floor(n/2);     
		r = n - 2*q;     
		B(i) = r; 
	end 
	B(i + 1) = q; 
	B = fliplr(B);
end