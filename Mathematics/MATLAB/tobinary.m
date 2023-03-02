function B = tobinary(n)

	i = 1; 	
	q = floor(n/2); 
	r = n - 2*q;
	B(i) = num2str(r(i)); 
	while 2 <= q     
		n = q;     
		i = i + 1;     
		q = floor(n/2);     
		r = n - 2*q;     
		B(i) = num2str(r); 
	end 
	B(i + 1) = num2str(q); 
	B = fliplr(B);
end