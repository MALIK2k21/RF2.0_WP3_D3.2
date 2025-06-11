function Total_grid_price = case2_GRID_PRICE(total_power)


Yearly_total_energy_consumption=total_power*1;

fixed_energy_price=2.86;  %sek/kWh
fixed_energy_price_euro=0.15; %Euro/kWh at 1 sek = 0.092 Euro at 5/27/2025



Total_grid_price=fixed_energy_price_euro*Yearly_total_energy_consumption;
end
